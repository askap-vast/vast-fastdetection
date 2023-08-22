FROM python:3.8

ENV \
  # python:
  PYTHONFAULTHANDLER=1 \
  PYTHONUNBUFFERED=1 \
  PYTHONHASHSEED=random \
  PYTHONDONTWRITEBYTECODE=1 \
  # pip:
  PIP_NO_CACHE_DIR=off \
  PIP_DISABLE_PIP_VERSION_CHECK=on \
  PIP_DEFAULT_TIMEOUT=100 \
  # poetry:
  POETRY_VERSION=1.1.12 \
  GET_POETRY_IGNORE_DEPRECATION=1 \
  POETRY_NO_INTERACTION=1 \
  POETRY_VIRTUALENVS_CREATE=false \
  POETRY_CACHE_DIR='/var/cache/pypoetry' \
  PATH="/root/.local/bin:$PATH"

# System dependencies:
RUN apt-get update && apt-get upgrade -y \
  && apt-get install --no-install-recommends -y \
    bash \
    build-essential \
    curl \
    gettext \
    git \
    wget \
    python3.8-dev \
  # Installing `poetry` package manager:
  && curl -sSL https://install.python-poetry.org | python3 \
  && poetry --version \
  # Cleaning cache:
  && apt-get purge -y --auto-remove -o APT::AutoRemove::RecommendsImportant=false \
  && apt-get clean -y && rm -rf /var/lib/apt/lists/*

WORKDIR /

COPY pyproject.toml poetry.lock /

# Install project dependencies
RUN poetry install --no-interaction --no-ansi

COPY . /

# Set the default command to run
#CMD ["poetry", "run", "python", "process_beam.py", "00"]
CMD ["python", "process_beam.py", "00"]
