from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json

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
    spin_enabled: bool
    softening: float = 0.055
    snapshot_times: tuple[float, ...] = (0.0, 2.25, 4.5)


def initialize_particles(
    cfg: Scenario,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(cfg.seed)
    radii = cfg.disk_radius * np.sqrt(rng.random(cfg.n))
    angles = 2.0 * np.pi * rng.random(cfg.n)
    pos = np.column_stack((radii * np.cos(angles), radii * np.sin(angles)))

    masses = np.full(cfg.n, cfg.total_mass / cfg.n)
    particle_radii = np.full(cfg.n, cfg.particle_radius)
    inertia = 0.5 * masses * particle_radii**2
    omega = np.zeros(cfg.n)

    r_safe = np.maximum(radii, 0.18 * cfg.disk_radius)
    tangential = np.column_stack((-np.sin(angles), np.cos(angles)))
    speeds = cfg.speed_scale / np.sqrt(r_safe)
    vel = tangential * speeds[:, None]
    vel += rng.normal(scale=0.025, size=vel.shape)

    vel -= np.average(vel, axis=0, weights=masses)
    pos -= np.average(pos, axis=0, weights=masses)

    return pos, vel, omega, masses, particle_radii, inertia


def compute_dynamics(
    pos: np.ndarray,
    vel: np.ndarray,
    omega: np.ndarray,
    masses: np.ndarray,
    radii: np.ndarray,
    inertia: np.ndarray,
    cfg: Scenario,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    n = len(masses)
    acc = np.zeros_like(pos)
    alpha = np.zeros_like(omega)
    grav_potential = 0.0
    rep_potential = 0.0

    for i in range(n - 1):
        rij = pos[i] - pos[i + 1 :]
        vij = vel[i] - vel[i + 1 :]
        dist2 = np.sum(rij * rij, axis=1) + cfg.softening**2
        dist = np.sqrt(dist2)
        unit = rij / dist[:, None]

        mass_prod = masses[i] * masses[i + 1 :]
        grav_force = (-mass_prod / (dist2 * dist))[:, None] * rij
        acc[i] += np.sum(grav_force, axis=0) / masses[i]
        acc[i + 1 :] -= grav_force / masses[i + 1 :, None]
        grav_potential += float(np.sum(-mass_prod / dist))

        overlap_limit = radii[i] + radii[i + 1 :]
        touching = dist < overlap_limit
        if not np.any(touching):
            continue

        contact_idx = np.nonzero(touching)[0]
        local_dist = dist[touching]
        local_unit = unit[touching]
        local_rel_vel = vij[touching]
        overlap = overlap_limit[touching]

        rep_mag = cfg.repulsion_k * ((overlap / local_dist) ** 8 - 1.0)
        rep_force = rep_mag[:, None] * local_unit

        tangent = np.column_stack((-local_unit[:, 1], local_unit[:, 0]))
        tangential_speed = np.sum(local_rel_vel * tangent, axis=1)

        if cfg.spin_enabled:
            slip_speed = tangential_speed - omega[i] * radii[i] - omega[i + 1 + contact_idx] * radii[i + 1 + contact_idx]
        else:
            slip_speed = tangential_speed

        friction_scalar = -cfg.friction_beta * slip_speed * rep_mag
        friction_force = friction_scalar[:, None] * tangent
        total_force = rep_force + friction_force

        acc[i] += np.sum(total_force, axis=0) / masses[i]
        acc[i + 1 + contact_idx] -= total_force / masses[i + 1 + contact_idx, None]

        if cfg.spin_enabled:
            alpha[i] += float(np.sum(-radii[i] * friction_scalar / inertia[i]))
            alpha[i + 1 + contact_idx] += -radii[i + 1 + contact_idx] * friction_scalar / inertia[i + 1 + contact_idx]

        rep_potential += float(
            np.sum(
                cfg.repulsion_k
                * ((overlap**8) / (7.0 * local_dist**7) + local_dist - 8.0 * overlap / 7.0)
            )
        )

    return acc, alpha, grav_potential, rep_potential


def translational_energy(vel: np.ndarray, masses: np.ndarray) -> float:
    return float(0.5 * np.sum(masses[:, None] * vel * vel))


def rotational_energy(omega: np.ndarray, inertia: np.ndarray) -> float:
    return float(0.5 * np.sum(inertia * omega * omega))


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


def angular_momentum(pos: np.ndarray, vel: np.ndarray, masses: np.ndarray, inertia: np.ndarray, omega: np.ndarray) -> tuple[float, float, float]:
    orbital = float(np.sum(masses * (pos[:, 0] * vel[:, 1] - pos[:, 1] * vel[:, 0])))
    spin = float(np.sum(inertia * omega))
    return orbital, spin, orbital + spin


def simulate(cfg: Scenario) -> dict[str, object]:
    pos, vel, omega, masses, radii, inertia = initialize_particles(cfg)
    acc, alpha, grav_pot, rep_pot = compute_dynamics(pos, vel, omega, masses, radii, inertia, cfg)

    times = np.arange(cfg.steps + 1) * cfg.dt
    energies = {
        "kinetic_trans": np.zeros(cfg.steps + 1),
        "kinetic_rot": np.zeros(cfg.steps + 1),
        "potential": np.zeros(cfg.steps + 1),
        "total": np.zeros(cfg.steps + 1),
    }
    metrics = {
        "rms_radius": np.zeros(cfg.steps + 1),
        "largest_cluster": np.zeros(cfg.steps + 1, dtype=int),
        "mean_abs_omega": np.zeros(cfg.steps + 1),
        "orbital_ang_momentum": np.zeros(cfg.steps + 1),
        "spin_ang_momentum": np.zeros(cfg.steps + 1),
        "total_ang_momentum": np.zeros(cfg.steps + 1),
    }
    snapshot_steps = sorted({min(cfg.steps, max(0, int(t / cfg.dt))) for t in cfg.snapshot_times})
    snapshots: dict[int, np.ndarray] = {}

    for step in range(cfg.steps + 1):
        pot = grav_pot + rep_pot
        kin_trans = translational_energy(vel, masses)
        kin_rot = rotational_energy(omega, inertia)
        orbital_l, spin_l, total_l = angular_momentum(pos, vel, masses, inertia, omega)

        energies["kinetic_trans"][step] = kin_trans
        energies["kinetic_rot"][step] = kin_rot
        energies["potential"][step] = pot
        energies["total"][step] = kin_trans + kin_rot + pot
        metrics["rms_radius"][step] = rms_radius(pos)
        metrics["largest_cluster"][step] = largest_cluster(pos, 3.2 * cfg.particle_radius)
        metrics["mean_abs_omega"][step] = float(np.mean(np.abs(omega)))
        metrics["orbital_ang_momentum"][step] = orbital_l
        metrics["spin_ang_momentum"][step] = spin_l
        metrics["total_ang_momentum"][step] = total_l
        if step in snapshot_steps:
            snapshots[step] = pos.copy()

        if step == cfg.steps:
            break

        pos_next = pos + vel * cfg.dt + 0.5 * acc * cfg.dt**2
        vel_half = vel + 0.5 * acc * cfg.dt
        omega_half = omega + 0.5 * alpha * cfg.dt
        acc_next, alpha_next, grav_pot, rep_pot = compute_dynamics(
            pos_next,
            vel_half,
            omega_half,
            masses,
            radii,
            inertia,
            cfg,
        )
        vel_next = vel + 0.5 * (acc + acc_next) * cfg.dt
        omega_next = omega + 0.5 * (alpha + alpha_next) * cfg.dt

        vel_next -= np.average(vel_next, axis=0, weights=masses)
        pos_next -= np.average(pos_next, axis=0, weights=masses)

        pos, vel, omega = pos_next, vel_next, omega_next
        acc, alpha = acc_next, alpha_next

    return {
        "config": cfg,
        "time": times,
        "energies": energies,
        "metrics": metrics,
        "snapshots": snapshots,
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


def moving_average(values: np.ndarray, window: int = 120) -> np.ndarray:
    kernel = np.ones(window) / window
    padded = np.pad(values, (window // 2, window - 1 - window // 2), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def plot_snapshots(baseline: dict[str, object], spin: dict[str, object]) -> None:
    cfg: Scenario = baseline["config"]  # type: ignore[assignment]
    steps = sorted(baseline["snapshots"].keys())  # type: ignore[assignment]
    fig, axes = plt.subplots(2, len(steps), figsize=(12.8, 7.4))
    titles = ["Старт", "Средний этап", "Финал"]

    for col, step in enumerate(steps):
        ax = axes[0, col]
        style_axes(ax, 1.45 * cfg.disk_radius)
        pts = baseline["snapshots"][step]  # type: ignore[index]
        ax.scatter(pts[:, 0], pts[:, 1], s=24, color="#8fd3ff", edgecolors="none", alpha=0.9)
        ax.set_title(f"{titles[col]}\nбез вращения", color="#222222", fontsize=12)

        ax = axes[1, col]
        style_axes(ax, 1.45 * cfg.disk_radius)
        pts = spin["snapshots"][step]  # type: ignore[index]
        ax.scatter(pts[:, 0], pts[:, 1], s=24, color="#ffb270", edgecolors="none", alpha=0.9)
        ax.set_title(f"{titles[col]}\nс вращением", color="#222222", fontsize=12)

    fig.tight_layout()
    fig.savefig(IMAGE_DIR / "spin_snapshots.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_energy_compare(baseline: dict[str, object], spin: dict[str, object]) -> None:
    time = baseline["time"]  # type: ignore[assignment]
    baseline_energy = baseline["energies"]  # type: ignore[assignment]
    spin_energy = spin["energies"]  # type: ignore[assignment]

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.7))

    axes[0].plot(time, baseline_energy["total"], label="Без вращения", color="#1f77b4", linewidth=2.1)
    axes[0].plot(time, spin_energy["total"], label="С вращением", color="#ff7f0e", linewidth=2.1)
    axes[0].set_title("Полная механическая энергия")
    axes[0].set_xlabel("Время")
    axes[0].set_ylabel("Энергия")
    axes[0].grid(alpha=0.3)
    axes[0].legend(frameon=False)

    axes[1].plot(time, spin_energy["kinetic_trans"], label="Поступательная", color="#1f77b4", linewidth=1.9)
    axes[1].plot(time, spin_energy["potential"], label="Потенциальная", color="#d62728", linewidth=1.9)
    axes[1].plot(time, spin_energy["total"], label="Полная", color="#111111", linewidth=2.2)
    axes[1].set_title("Разложение энергии при учёте вращения")
    axes[1].set_xlabel("Время")
    axes[1].set_ylabel("Энергия")
    axes[1].grid(alpha=0.3)
    rot_ax = axes[1].twinx()
    rot_ax.plot(time, spin_energy["kinetic_rot"], label="Вращательная", color="#2ca02c", linewidth=1.9)
    rot_ax.set_ylabel("Вращательная энергия", color="#2ca02c")
    rot_ax.tick_params(axis="y", colors="#2ca02c")

    lines = axes[1].get_lines() + rot_ax.get_lines()
    labels = [line.get_label() for line in lines]
    axes[1].legend(lines, labels, frameon=False, loc="lower left")

    fig.tight_layout()
    fig.savefig(IMAGE_DIR / "spin_energy_compare.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_metrics(baseline: dict[str, object], spin: dict[str, object]) -> None:
    time = baseline["time"]  # type: ignore[assignment]
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.7))

    axes[0].plot(time, baseline["metrics"]["rms_radius"], label="Без вращения", color="#1f77b4", linewidth=2.0)  # type: ignore[index]
    axes[0].plot(time, spin["metrics"]["rms_radius"], label="С вращением", color="#ff7f0e", linewidth=2.0)  # type: ignore[index]
    axes[0].set_title("Среднеквадратичный радиус диска")
    axes[0].set_xlabel("Время")
    axes[0].set_ylabel(r"$R_{\mathrm{rms}}$")
    axes[0].grid(alpha=0.3)
    axes[0].legend(frameon=False)

    baseline_cluster = moving_average(baseline["metrics"]["largest_cluster"].astype(float))  # type: ignore[index]
    spin_cluster = moving_average(spin["metrics"]["largest_cluster"].astype(float))  # type: ignore[index]
    axes[1].plot(time, baseline_cluster, label="Без вращения", color="#1f77b4", linewidth=2.0)
    axes[1].plot(time, spin_cluster, label="С вращением", color="#ff7f0e", linewidth=2.0)
    axes[1].set_title("Размер крупнейшего кластера")
    axes[1].set_xlabel("Время")
    axes[1].set_ylabel("Частиц")
    axes[1].grid(alpha=0.3)
    axes[1].legend(frameon=False)

    fig.tight_layout()
    fig.savefig(IMAGE_DIR / "spin_metrics.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_angular_momentum(spin: dict[str, object]) -> None:
    time = spin["time"]  # type: ignore[assignment]
    metrics = spin["metrics"]  # type: ignore[assignment]
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.7))

    axes[0].plot(time, metrics["orbital_ang_momentum"], label="Орбитальный", color="#1f77b4", linewidth=1.8)
    axes[0].plot(time, metrics["total_ang_momentum"], label="Суммарный", color="#111111", linewidth=2.2)
    axes[0].set_title("Момент импульса")
    axes[0].set_xlabel("Время")
    axes[0].set_ylabel(r"$L_z$")
    axes[0].grid(alpha=0.3)
    spin_ax = axes[0].twinx()
    spin_ax.plot(time, metrics["spin_ang_momentum"], label="Спиновый", color="#2ca02c", linewidth=1.8)
    spin_ax.set_ylabel(r"$L^{\mathrm{spin}}_z$", color="#2ca02c")
    spin_ax.tick_params(axis="y", colors="#2ca02c")

    lines = axes[0].get_lines() + spin_ax.get_lines()
    labels = [line.get_label() for line in lines]
    axes[0].legend(lines, labels, frameon=False, loc="center right")

    axes[1].plot(time, metrics["mean_abs_omega"], color="#9467bd", linewidth=2.0)
    axes[1].set_title("Средняя по модулю угловая скорость")
    axes[1].set_xlabel("Время")
    axes[1].set_ylabel(r"$\langle |\omega| \rangle$")
    axes[1].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(IMAGE_DIR / "spin_angular_momentum.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_summary(baseline: dict[str, object], spin: dict[str, object]) -> None:
    def summarize(sim: dict[str, object]) -> dict[str, float]:
        energies = sim["energies"]  # type: ignore[assignment]
        metrics = sim["metrics"]  # type: ignore[assignment]
        total = energies["total"]
        total_l = metrics["total_ang_momentum"]
        rot = energies["kinetic_rot"]
        return {
            "energy_start": float(total[0]),
            "energy_end": float(total[-1]),
            "energy_relative_change": float((total[-1] - total[0]) / max(abs(total[0]), 1e-9)),
            "rms_start": float(metrics["rms_radius"][0]),
            "rms_end": float(metrics["rms_radius"][-1]),
            "cluster_start": int(metrics["largest_cluster"][0]),
            "cluster_end": int(metrics["largest_cluster"][-1]),
            "mean_abs_omega_end": float(metrics["mean_abs_omega"][-1]),
            "rotational_energy_end": float(rot[-1]),
            "angular_momentum_relative_change": float((total_l[-1] - total_l[0]) / max(abs(total_l[0]), 1e-9)),
        }

    summary = {
        "baseline_no_spin": summarize(baseline),
        "spin_enabled": summarize(spin),
    }
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


def main() -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)

    base_cfg = dict(
        n=96,
        steps=1800,
        dt=0.0025,
        seed=314,
        disk_radius=1.0,
        total_mass=0.5,
        particle_radius=0.042,
        speed_scale=0.35,
        repulsion_k=0.0018,
        friction_beta=0.7,
        softening=0.055,
        snapshot_times=(0.0, 2.25, 4.5),
    )
    baseline_cfg = Scenario(
        name="baseline_no_spin",
        label="Самогравитация с трением без вращения",
        spin_enabled=False,
        **base_cfg,
    )
    spin_cfg = Scenario(
        name="spin_enabled",
        label="Самогравитация с трением и собственным вращением",
        spin_enabled=True,
        **base_cfg,
    )

    baseline = simulate(baseline_cfg)
    spin = simulate(spin_cfg)

    plot_snapshots(baseline, spin)
    plot_energy_compare(baseline, spin)
    plot_metrics(baseline, spin)
    plot_angular_momentum(spin)
    write_summary(baseline, spin)


if __name__ == "__main__":
    main()
