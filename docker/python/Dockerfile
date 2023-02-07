FROM python:3.8

RUN pip install --upgrade pip

RUN mkdir /app
COPY . /app
#COPY pyproject.toml /app

WORKDIR /app

ENV PYTHONPATH=${PYTHONPATH}:${PWD}

RUN #pip3 --upgrade pip
RUN pip3 install pyregion
RUN pip3 install pyregion

#RUN poetry config virtualenvs.create false
#RUN poetry add aplpy
