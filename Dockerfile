FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONFAULTHANDLER=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONHASHSEED=random \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=off \
    PIP_DISABLE_PIP_VERSION_CHECK=on \
    PIP_DEFAULT_TIMEOUT=100 \
    POETRY_VERSION=1.1.12 \
    POETRY_NO_INTERACTION=1 \
    POETRY_CACHE_DIR='/var/cache/pypoetry' \
    PATH="/root/.local/bin:$PATH"

# --- Python 3.10 from deadsnakes ---
RUN apt-get update && apt-get install -y --no-install-recommends \
      software-properties-common ca-certificates gnupg \
 && add-apt-repository -y ppa:deadsnakes/ppa \
 && apt-get update && apt-get install -y --no-install-recommends \
      python3.10 python3.10-dev python3.10-venv python3.10-distutils \
 && rm -rf /var/lib/apt/lists/*

# --- System build deps ---
RUN apt-get update && apt-get install -y --no-install-recommends \
      bash build-essential curl gettext git wget gfortran pkg-config \
 && rm -rf /var/lib/apt/lists/*

# --- Pip tooling for Python 3.10 ---
RUN python3.10 -m ensurepip --upgrade \
 && python3.10 -m pip install --upgrade \
      "pip==23.2.1" \
      "wheel<0.41" \
      "setuptools==58.5.3" \
      "packaging>=24.1"

# --- Poetry under Python 3.10 ---
RUN python3.10 -m pip install "poetry==${POETRY_VERSION}" \
 && poetry --version

WORKDIR /code/

# --- Copy dependency files first ---
COPY pyproject.toml poetry.lock /code/

# --- Install project deps inside a venv ---
ENV PIP_NO_BUILD_ISOLATION=1
RUN poetry config virtualenvs.create true \
 && poetry env use /usr/bin/python3.10 \
 && poetry run python -V \
 && poetry run python -m pip install -U "pip==23.2.1" "wheel<0.41" "setuptools==58.5.3" \
 && poetry install --no-interaction --no-ansi \
 # expose venv globally so "python" works
 && VENV_DIR="$(poetry env info --path)" \
 && ln -s "$VENV_DIR" /opt/venv
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# --- Copy the rest of the source ---
COPY . /code/

# --- Default command ---
CMD ["python", "/code/process_beam.py", "00"]

