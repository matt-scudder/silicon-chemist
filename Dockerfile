FROM python:3.13-slim-bookworm

WORKDIR /usr/src/app

RUN apt update && apt install -y --no-install-recommends \
    default-jre \
    swig \
    libopenbabel-dev \
    g++ \
    && rm -rf /var/lib/apt/lists/*
RUN ln -s /usr/include/openbabel3 /usr/local/include/openbabel3

COPY requirements.txt .
RUN --mount=type=cache,target=/root/.cache/pip pip install -r requirements.txt

COPY . .

EXPOSE 5000

CMD ["python3", "runserver.py"]