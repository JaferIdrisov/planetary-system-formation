#!/bin/bash

# Скрипт инициализации Git репозитория для проекта
# Запустите этот скрипт после создания репозитория на GitHub

echo "=== Инициализация Git репозитория ==="

cd "$(dirname "$0")"

# Инициализация Git
git init
git branch -M main

# Добавление всех файлов
git add .
git commit -m "Initial commit: сайт проекта Образование планетной системы"

echo ""
echo "=== Следующие шаги ==="
echo "1. Создайте репозиторий на GitHub с именем 'planetary-system-formation'"
echo "2. Выполните команды:"
echo "   git remote add origin https://github.com/daidrisov/planetary-system-formation.git"
echo "   git push -u origin main"
echo ""
echo "3. В настройках репозитория на GitHub:"
echo "   - Перейдите в Settings > Pages"
echo "   - В разделе 'Build and deployment' выберите 'GitHub Actions'"
echo ""
echo "4. После пуша workflow автоматически запустится и развернёт сайт"
