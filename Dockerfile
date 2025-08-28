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

# --- Python 3.10 (deadsnakes) ---
RUN apt-get update && apt-get install -y --no-install-recommends \
      software-properties-common ca-certificates gnupg \
 && add-apt-repository -y ppa:deadsnakes/ppa \
 && apt-get update && apt-get install -y --no-install-recommends \
      python3.10 python3.10-dev python3.10-venv python3.10-distutils \
 && rm -rf /var/lib/apt/lists/*

# Base dev tools (similar to your original)
RUN apt-get update && apt-get install -y --no-install-recommends \
      bash build-essential curl gettext git wget \
 && rm -rf /var/lib/apt/lists/*

# Pip/tooling for py3.10 (use known-good versions)
RUN python3.10 -m ensurepip --upgrade \
 && python3.10 -m pip install --upgrade \
      "pip==23.2.1" \
      "wheel<0.41" \
      "setuptools<75" \
      "packaging>=24.1"

# Install Poetry for py3.10 (no curl-installer)
RUN python3.10 -m pip install "poetry==${POETRY_VERSION}" \
 && poetry --version

WORKDIR /code/

# Copy dependency files first (cache-friendly)
COPY pyproject.toml poetry.lock /code/

# IMPORTANT: avoid pip build isolation so legacy packages (pyregion) use our pinned setuptools
ENV PIP_NO_BUILD_ISOLATION=1

# Use a Poetry-managed venv (Python 3.10), pre-pin build tools *inside* that venv, then install deps
RUN poetry config virtualenvs.create true \
 && poetry env use /usr/bin/python3.10 \
 && poetry run python -V \
 # old setuptools API needed by pyregion 2.1.1 (and some astropy-helpers builds)
 && poetry run python -m pip install -U "pip==23.2.1" "wheel<0.41" "setuptools==58.5.3" \
 # (optional) if pyregion needs cython on your stack, uncomment the next line
 # && poetry run python -m pip install "cython<3" \
 && poetry install --no-interaction --no-ansi

# Bring in the rest of the source
COPY . /code/

# Default command
CMD ["poetry", "run", "python", "/code/process_beam.py", "00"]

