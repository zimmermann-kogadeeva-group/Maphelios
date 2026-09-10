# app/Dockerfile

FROM python:3.12-slim

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

ENV STREAMLIT_SERVER_MAX_UPLOAD_SIZE=1000

WORKDIR /app

COPY pyproject.toml /app
COPY app/ /app/app/
COPY src/ /app/src/

RUN uv sync --extra streamlit --no-dev --no-cache

# Installing deps
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    zip

# Adding BLAST
RUN curl https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-x64-linux.tar.gz -o blast.tar.gz && \
    tar -C /opt/ -zxvf blast.tar.gz && \
    rm blast.tar.gz

RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.31/minimap2-2.31_x64-linux.tar.bz2 -o minimap2.tar.bz2 && \
    tar -C /opt/ -jxvf minimap2.tar.bz2 && \
    rm minimap2.tar.bz2

# Updating path env var
ENV PATH="/app/.venv/bin:/opt/ncbi-blast-2.16.0+/bin:/opt/minimap2-2.31_x64_linux:${PATH}"

# Cleaning up
RUN apt remove --purge -y curl zip && apt clean && rm -rf /var/lib/apt/lists/*

EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["streamlit", "run", "app/app.py", "--server.port=8501", "--server.address=0.0.0.0"]
