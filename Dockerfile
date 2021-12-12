FROM ubuntu:20.04
LABEL maintainer="https://github.com/pdiakumis"

ARG MINI_VERSION=4.11.0-0
ARG MINI_URL=https://github.com/conda-forge/miniforge/releases/download/${MINI_VERSION}/Mambaforge-${MINI_VERSION}-Linux-x86_64.sh

# install core pkgs, mambaforge
RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
    bash bzip2 curl git less vim wget zip ca-certificates && \
    apt-get clean && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    curl -L "${MINI_URL}" -o "mambaforge.sh" && \
    /bin/bash mambaforge.sh -b -p /opt/mambaforge/ && \
    rm mambaforge.sh

# create conda envs
ENV PATH="/opt/mambaforge/bin:$PATH"
ARG PCGR_CONDA_ENV_DIR=/home/pcgr_conda_envs
COPY ./conda/env ${PCGR_CONDA_ENV_DIR}
RUN mamba env create --file ${PCGR_CONDA_ENV_DIR}/pcgr.yml
RUN mamba env create --file ${PCGR_CONDA_ENV_DIR}/pcgrr.yml
RUN mamba clean --all --force-pkgs-dirs --yes

ARG PCGR_ENV_NAME="pcgr"
# pcgr env is activated by default
ENV PATH="/opt/mambaforge/envs/${PCGR_ENV_NAME}/bin:${PATH}"
ENV CONDA_PREFIX="/opt/mambaforge/envs/${PCGR_ENV_NAME}"
