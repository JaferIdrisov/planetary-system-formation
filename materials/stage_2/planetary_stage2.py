from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import math

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
IMAGE_DIR = ROOT / "image"
SUMMARY_PATH = ROOT / "summary.json"


@dataclass(frozen=True)
class Scenario:
    name: str
    label: str
    n: int
    steps: int
    dt: float
    seed: int
    disk_radius: float
    total_mass: float
    particle_radius: float
    speed_scale: float
    repulsion_k: float
    friction_beta: float
    central_mass: float = 0.0
    pairwise_gravity: bool = True
    softening: float = 0.025
    snapshot_times: tuple[float, ...] = (0.0, 2.0, 4.0)


def initialize_particles(cfg: Scenario) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(cfg.seed)
    radii = cfg.disk_radius * np.sqrt(rng.random(cfg.n))
    angles = 2.0 * np.pi * rng.random(cfg.n)
    pos = np.column_stack((radii * np.cos(angles), radii * np.sin(angles)))
    masses = np.full(cfg.n, cfg.total_mass / cfg.n)
    particle_radii = np.full(cfg.n, cfg.particle_radius)

    r_safe = np.maximum(radii, 0.18 * cfg.disk_radius)
    tangential = np.column_stack((-np.sin(angles), np.cos(angles)))
    speeds = cfg.speed_scale / np.sqrt(r_safe)
    vel = tangential * speeds[:, None]
    vel += rng.normal(scale=0.025, size=vel.shape)

    if cfg.central_mass == 0.0:
        vel -= np.average(vel, axis=0, weights=masses)
        pos -= np.average(pos, axis=0, weights=masses)

    return pos, vel, masses, particle_radii


def compute_forces(
    pos: np.ndarray,
    vel: np.ndarray,
    masses: np.ndarray,
    radii: np.ndarray,
    cfg: Scenario,
) -> tuple[np.ndarray, float, float]:
    n = len(masses)
    acc = np.zeros_like(pos)
    grav_potential = 0.0
    rep_potential = 0.0

    if cfg.central_mass > 0.0:
        r2 = np.sum(pos * pos, axis=1) + cfg.softening**2
        r = np.sqrt(r2)
        acc += (-cfg.central_mass * pos.T / (r2 * r)).T
        grav_potential += float(np.sum(-cfg.central_mass * masses / r))

    if not cfg.pairwise_gravity and cfg.repulsion_k <= 0.0:
        return acc, grav_potential, rep_potential

    for i in range(n - 1):
        rij = pos[i] - pos[i + 1 :]
        vij = vel[i] - vel[i + 1 :]
        dist2 = np.sum(rij * rij, axis=1) + cfg.softening**2
        dist = np.sqrt(dist2)
        unit = rij / dist[:, None]

        if cfg.pairwise_gravity:
            mass_prod = masses[i] * masses[i + 1 :]
            grav_force = (-mass_prod / (dist2 * dist))[:, None] * rij
            acc[i] += np.sum(grav_force, axis=0) / masses[i]
            acc[i + 1 :] -= grav_force / masses[i + 1 :, None]
            grav_potential += float(np.sum(-mass_prod / dist))

        if cfg.repulsion_k <= 0.0:
            continue

        overlap_limit = radii[i] + radii[i + 1 :]
        touching = dist < overlap_limit
        if not np.any(touching):
            continue

        contact_idx = np.nonzero(touching)[0]
        overlap = overlap_limit[touching]
        local_dist = dist[touching]
        local_unit = unit[touching]
        local_rel_vel = vij[touching]

        rep_mag = cfg.repulsion_k * ((overlap / local_dist) ** 8 - 1.0)
        rep_force = rep_mag[:, None] * local_unit

        tangent = np.column_stack((-local_unit[:, 1], local_unit[:, 0]))
        tangential_speed = np.sum(local_rel_vel * tangent, axis=1)
        friction_force = (-cfg.friction_beta * tangential_speed * rep_mag)[:, None] * tangent
        total_force = rep_force + friction_force

        acc[i] += np.sum(total_force, axis=0) / masses[i]
        acc[i + 1 + contact_idx] -= total_force / masses[i + 1 + contact_idx, None]
        rep_potential += float(
            np.sum(
                cfg.repulsion_k
                * ((overlap**8) / (7.0 * local_dist**7) + local_dist - 8.0 * overlap / 7.0)
            )
        )

    return acc, grav_potential, rep_potential


def kinetic_energy(vel: np.ndarray, masses: np.ndarray) -> float:
    return float(0.5 * np.sum(masses[:, None] * vel * vel))


def rms_radius(pos: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum(pos * pos, axis=1))))


def largest_cluster(pos: np.ndarray, threshold: float) -> int:
    n = len(pos)
    dist2 = np.sum((pos[:, None, :] - pos[None, :, :]) ** 2, axis=2)
    adjacency = dist2 < threshold**2
    visited = np.zeros(n, dtype=bool)
    best = 0

    for start in range(n):
        if visited[start]:
            continue
        stack = [start]
        visited[start] = True
        size = 0
        while stack:
            node = stack.pop()
            size += 1
            neighbors = np.nonzero(adjacency[node])[0]
            for nxt in neighbors:
                if not visited[nxt]:
                    visited[nxt] = True
                    stack.append(int(nxt))
        best = max(best, size)
    return best


def simulate(cfg: Scenario) -> dict[str, object]:
    pos, vel, masses, radii = initialize_particles(cfg)
    acc, grav_pot, rep_pot = compute_forces(pos, vel, masses, radii, cfg)

    times = np.arange(cfg.steps + 1) * cfg.dt
    energies = {
        "kinetic": np.zeros(cfg.steps + 1),
        "potential": np.zeros(cfg.steps + 1),
        "total": np.zeros(cfg.steps + 1),
    }
    metrics = {
        "rms_radius": np.zeros(cfg.steps + 1),
        "largest_cluster": np.zeros(cfg.steps + 1, dtype=int),
    }
    track_count = min(18, cfg.n)
    tracks = np.zeros((cfg.steps + 1, track_count, 2))
    snapshot_steps = sorted({min(cfg.steps, max(0, int(t / cfg.dt))) for t in cfg.snapshot_times})
    snapshots: dict[int, np.ndarray] = {}

    for step in range(cfg.steps + 1):
        pot = grav_pot + rep_pot
        kin = kinetic_energy(vel, masses)
        energies["kinetic"][step] = kin
        energies["potential"][step] = pot
        energies["total"][step] = kin + pot
        metrics["rms_radius"][step] = rms_radius(pos)
        metrics["largest_cluster"][step] = largest_cluster(pos, 3.2 * cfg.particle_radius)
        tracks[step] = pos[:track_count]
        if step in snapshot_steps:
            snapshots[step] = pos.copy()

        if step == cfg.steps:
            break

        pos_next = pos + vel * cfg.dt + 0.5 * acc * cfg.dt**2
        vel_half = vel + 0.5 * acc * cfg.dt
        acc_next, grav_pot, rep_pot = compute_forces(pos_next, vel_half, masses, radii, cfg)
        vel_next = vel + 0.5 * (acc + acc_next) * cfg.dt

        if cfg.central_mass == 0.0:
            vel_next -= np.average(vel_next, axis=0, weights=masses)
            pos_next -= np.average(pos_next, axis=0, weights=masses)

        pos, vel, acc = pos_next, vel_next, acc_next

    return {
        "config": cfg,
        "time": times,
        "energies": energies,
        "metrics": metrics,
        "snapshots": snapshots,
        "tracks": tracks,
    }


def style_axes(ax: plt.Axes, limit: float) -> None:
    ax.set_aspect("equal")
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_facecolor("#07111f")
    ax.grid(color="#21344d", alpha=0.3, linewidth=0.5)
    ax.tick_params(colors="#d5dde5", labelsize=8)
    for spine in ax.spines.values():
        spine.set_color("#516377")


def plot_central_orbits(sim: dict[str, object]) -> None:
    cfg: Scenario = sim["config"]  # type: ignore[assignment]
    tracks = sim["tracks"]  # type: ignore[assignment]
    fig, ax = plt.subplots(figsize=(7.4, 7.0))
    style_axes(ax, 1.8 * cfg.disk_radius)
    colors = plt.cm.viridis(np.linspace(0.15, 0.95, tracks.shape[1]))
    for idx, color in enumerate(colors):
        ax.plot(tracks[:, idx, 0], tracks[:, idx, 1], color=color, linewidth=1.1, alpha=0.9)
        ax.scatter(tracks[-1, idx, 0], tracks[-1, idx, 1], color=color, s=22, edgecolors="none")
    ax.scatter(0.0, 0.0, s=160, color="#ffcf5c", edgecolors="white", linewidth=0.8, zorder=3)
    ax.set_title("Проверка сценария с центральной звездой", color="#222222", fontsize=13)
    ax.set_xlabel("x", color="#d5dde5")
    ax.set_ylabel("y", color="#d5dde5")
    fig.tight_layout()
    fig.savefig(IMAGE_DIR / "central_orbits.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_energy(sim: dict[str, object], filename: str, title: str) -> None:
    time = sim["time"]  # type: ignore[assignment]
    energies = sim["energies"]  # type: ignore[assignment]

    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    ax.plot(time, energies["kinetic"], label="Кинетическая", color="#1f77b4", linewidth=1.8)
    ax.plot(time, energies["potential"], label="Потенциальная", color="#d62728", linewidth=1.8)
    ax.plot(time, energies["total"], label="Полная", color="#111111", linewidth=2.2)
    ax.set_xlabel("Время")
    ax.set_ylabel("Энергия")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(IMAGE_DIR / filename, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_disk_snapshots(no_friction: dict[str, object], with_friction: dict[str, object]) -> None:
    cfg: Scenario = no_friction["config"]  # type: ignore[assignment]
    steps = sorted(no_friction["snapshots"].keys())  # type: ignore[assignment]
    fig, axes = plt.subplots(2, len(steps), figsize=(12.8, 7.6))
    titles = ["Старт", "Средний этап", "Финал"]

    for col, step in enumerate(steps):
        ax = axes[0, col]
        style_axes(ax, 1.45 * cfg.disk_radius)
        pts = no_friction["snapshots"][step]  # type: ignore[index]
        ax.scatter(pts[:, 0], pts[:, 1], s=24, color="#8fd3ff", edgecolors="none", alpha=0.9)
        ax.set_title(f"{titles[col]}\nбез трения", color="#222222", fontsize=12)

        ax = axes[1, col]
        style_axes(ax, 1.45 * cfg.disk_radius)
        pts = with_friction["snapshots"][step]  # type: ignore[index]
        ax.scatter(pts[:, 0], pts[:, 1], s=24, color="#ffb270", edgecolors="none", alpha=0.9)
        ax.set_title(f"{titles[col]}\nс трением", color="#222222", fontsize=12)

    fig.tight_layout()
    fig.savefig(IMAGE_DIR / "disk_snapshots.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_aggregation_metrics(no_friction: dict[str, object], with_friction: dict[str, object]) -> None:
    time = no_friction["time"]  # type: ignore[assignment]
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.6))

    def moving_average(values: np.ndarray, window: int = 120) -> np.ndarray:
        kernel = np.ones(window) / window
        padded = np.pad(values, (window // 2, window - 1 - window // 2), mode="edge")
        return np.convolve(padded, kernel, mode="valid")

    no_friction_cluster = moving_average(no_friction["metrics"]["largest_cluster"].astype(float))  # type: ignore[index]
    with_friction_cluster = moving_average(with_friction["metrics"]["largest_cluster"].astype(float))  # type: ignore[index]

    axes[0].plot(time, no_friction["metrics"]["rms_radius"], label="Без трения", color="#1f77b4", linewidth=2.0)  # type: ignore[index]
    axes[0].plot(time, with_friction["metrics"]["rms_radius"], label="С трением", color="#ff7f0e", linewidth=2.0)  # type: ignore[index]
    axes[0].set_title("Среднеквадратичный радиус")
    axes[0].set_xlabel("Время")
    axes[0].set_ylabel(r"$R_{\mathrm{rms}}$")
    axes[0].grid(alpha=0.3)
    axes[0].legend(frameon=False)

    axes[1].plot(time, no_friction_cluster, label="Без трения", color="#1f77b4", linewidth=2.0)
    axes[1].plot(time, with_friction_cluster, label="С трением", color="#ff7f0e", linewidth=2.0)
    axes[1].set_title("Скользящее среднее размера крупнейшего кластера")
    axes[1].set_xlabel("Время")
    axes[1].set_ylabel("Частиц")
    axes[1].grid(alpha=0.3)
    axes[1].legend(frameon=False)
    axes[1].scatter(time[-1], no_friction_cluster[-1], color="#1f77b4", s=28, zorder=3)
    axes[1].scatter(time[-1], with_friction_cluster[-1], color="#ff7f0e", s=28, zorder=3)
    axes[1].annotate(
        f"финал: {int(no_friction['metrics']['largest_cluster'][-1])}",  # type: ignore[index]
        (time[-1], no_friction_cluster[-1]),
        xytext=(-82, -18),
        textcoords="offset points",
        color="#1f77b4",
        fontsize=9,
    )
    axes[1].annotate(
        f"финал: {int(with_friction['metrics']['largest_cluster'][-1])}",  # type: ignore[index]
        (time[-1], with_friction_cluster[-1]),
        xytext=(-82, 8),
        textcoords="offset points",
        color="#ff7f0e",
        fontsize=9,
    )

    fig.tight_layout()
    fig.savefig(IMAGE_DIR / "aggregation_metrics.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_summary(central: dict[str, object], no_friction: dict[str, object], with_friction: dict[str, object]) -> None:
    def summarize(sim: dict[str, object]) -> dict[str, float]:
        energies = sim["energies"]  # type: ignore[assignment]
        metrics = sim["metrics"]  # type: ignore[assignment]
        total = energies["total"]
        return {
            "energy_start": float(total[0]),
            "energy_end": float(total[-1]),
            "energy_relative_drift": float((total[-1] - total[0]) / max(abs(total[0]), 1e-9)),
            "rms_start": float(metrics["rms_radius"][0]),
            "rms_end": float(metrics["rms_radius"][-1]),
            "cluster_start": int(metrics["largest_cluster"][0]),
            "cluster_end": int(metrics["largest_cluster"][-1]),
        }

    summary = {
        "central": summarize(central),
        "disk_no_friction": summarize(no_friction),
        "disk_with_friction": summarize(with_friction),
    }
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


def main() -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)

    central_cfg = Scenario(
        name="central",
        label="Центральная звезда",
        n=48,
        steps=900,
        dt=0.01,
        seed=2026,
        disk_radius=1.45,
        total_mass=0.12,
        particle_radius=0.028,
        speed_scale=0.95,
        repulsion_k=0.0,
        friction_beta=0.0,
        central_mass=1.0,
        pairwise_gravity=False,
        snapshot_times=(0.0, 4.5, 9.0),
    )
    disk_base = dict(
        n=96,
        steps=1800,
        dt=0.0025,
        seed=314,
        disk_radius=1.0,
        total_mass=0.5,
        particle_radius=0.042,
        speed_scale=0.35,
        repulsion_k=0.0018,
        softening=0.055,
        snapshot_times=(0.0, 2.25, 4.5),
    )
    disk_no_friction_cfg = Scenario(
        name="disk_no_friction",
        label="Самогравитация без трения",
        friction_beta=0.0,
        **disk_base,
    )
    disk_with_friction_cfg = Scenario(
        name="disk_with_friction",
        label="Самогравитация с трением",
        friction_beta=0.7,
        **disk_base,
    )

    central = simulate(central_cfg)
    no_friction = simulate(disk_no_friction_cfg)
    with_friction = simulate(disk_with_friction_cfg)

    plot_central_orbits(central)
    plot_energy(central, "central_energy.png", "Энергия в сценарии с центральной звездой")
    plot_disk_snapshots(no_friction, with_friction)
    plot_energy(no_friction, "disk_energy_no_friction.png", "Энергия самогравитирующего диска без трения")
    plot_energy(with_friction, "disk_energy_with_friction.png", "Энергия самогравитирующего диска с трением")
    plot_aggregation_metrics(no_friction, with_friction)
    write_summary(central, no_friction, with_friction)


if __name__ == "__main__":
    main()
