# Сайт проекта "Образование планетной системы"

Сайт группового проекта, созданный на основе шаблона [HugoBlox theme-research-group](https://github.com/HugoBlox/theme-research-group).

## Рабочая группа

- **Идрисов Джафер Арсенович** - 1132232876
- **Ведьмина Александра Сергеевна** - 1132236003
- **Кузьмин Егор Витальевич** - 1132236046
- **Курилко-Рюмин Евгений Михайлович** - 1132232883
- **Закиров Нурислам Дамирович** - 1132236040
- **Глущенко Евгений Игоревич** - 1132239110

## Этапы проекта

1. **Модель** - Презентация по научной проблеме. Теоретическое описание задачи. Описание модели.
2. **Алгоритмы** - Презентация по алгоритмам решения задачи.
3. **Комплексы программ** - Описание программной реализации проекта.
4. **Защита проекта** - Коллективное обсуждение результата проекта, самооценка деятельности.

## Быстрый старт

### Вариант 1: Использование Hugo (рекомендуется)

1. **Установите Hugo Extended** (версия 0.100 или выше):
   
   **Ubuntu/Debian:**
   ```bash
   sudo apt-get install hugo
   ```
   
   **macOS:**
   ```bash
   brew install hugo
   ```
   
   **Windows:**
   ```bash
   choco install hugo-extended
   ```
   
   Или скачайте с [официального сайта](https://gohugo.io/getting-started/installing/)

2. **Перейдите в директорию сайта:**
   ```bash
   cd site
   ```

3. **Установите зависимости Go:**
   ```bash
   go mod init github.com/daidrisov/planetary-system-formation
   go mod tidy
   ```

4. **Запустите локальный сервер:**
   ```bash
   hugo server
   ```

5. **Откройте браузер** по адресу: http://localhost:1313

### Вариант 2: Docker

Если у вас установлен Docker:

```bash
cd site
docker run --rm -it -p 1313:1313 -v $(pwd):/src klakegg/hugo:ext-alpine server
```

## Развёртывание на GitHub Pages

### Шаг 1: Создайте репозиторий на GitHub

1. Создайте новый репозиторий с именем `planetary-system-formation`
2. Не инициализируйте его README, .gitignore или лицензией

### Шаг 2: Инициализируйте Git репозиторий

```bash
cd site
git init
git add .
git commit -m "Initial commit"
git branch -M main
git remote add origin https://github.com/daidrisov/planetary-system-formation.git
git push -u origin main
```

### Шаг 3: Настройте GitHub Pages

1. Перейдите в настройки репозитория на GitHub
2. Выберите раздел **Pages**
3. В разделе **Build and deployment**:
   - Source: GitHub Actions
4. Workflow автоматически запустится после пуша

### Шаг 4: Проверьте статус развёртывания

1. Перейдите во вкладку **Actions** вашего репозитория
2. Дождитесь завершения workflow
3. После успешного развёртывания сайт будет доступен по адресу:
   `https://daidrisov.github.io/planetary-system-formation/`

## Структура проекта

```
site/
├── config/              # Конфигурация Hugo
├── content/             # Контент сайта
│   ├── authors/         # Профили участников
│   ├── people/          # Страница команды
│   ├── stages/          # Этапы проекта
│   └── contact/         # Контакты
├── static/              # Статические файлы
│   └── uploads/         # Загруженные файлы (презентации, отчеты)
├── .github/workflows/   # GitHub Actions workflow
└── theme.toml           # Информация о теме
```

## Добавление материалов

### Добавление презентации или отчета

1. Поместите PDF файл в `static/uploads/`
2. Обновите соответствующую страницу в `content/stages/`

### Обновление информации об участниках

Отредактируйте файлы в `content/authors/`:
- `idrisov-dzhafar/_index.md`
- `vedmina-aleksandra/_index.md`
- `kuzmin-egor/_index.md`
- `kurilko-ryumin-evgeniy/_index.md`
- `zakirov-nurislam/_index.md`
- `glushhenko-evgeniy/_index.md`

## Лицензия

Этот проект использует тему [HugoBlox theme-research-group](https://github.com/HugoBlox/theme-research-group) под лицензией MIT.
