FROM ubuntu:20.04
LABEL org.opencontainers.image.authors='sigven@ifi.uio.no, peterdiakumis@gmail.com' \
      org.opencontainers.image.description='Personal Cancer Genome Reporter (PCGR)' \
      org.opencontainers.image.source='https://github.com/sigven/pcgr' \
      org.opencontainers.image.url='https://github.com/sigven/pcgr' \
      org.opencontainers.image.documentation='https://sigven.github.io/pcgr' \
      org.opencontainers.image.licenses='MIT'

ARG MINI_VERSION=4.11.0-0
ARG MINI_URL=https://github.com/conda-forge/miniforge/releases/download/${MINI_VERSION}/Mambaforge-${MINI_VERSION}-Linux-x86_64.sh

# install core pkgs, mambaforge
RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
    bash bzip2 curl less wget zip ca-certificates && \
    apt-get clean && \
    curl --silent -L "${MINI_URL}" -o "mambaforge.sh" && \
    /bin/bash mambaforge.sh -b -p /opt/mambaforge/ && \
    rm mambaforge.sh

# create conda envs
ENV PATH="/opt/mambaforge/bin:$PATH"
ARG PCGR_CONDA_ENV_DIR=/home/pcgr_conda_envs
COPY ./conda/env/lock ${PCGR_CONDA_ENV_DIR}
RUN mamba create -n pcgr --file ${PCGR_CONDA_ENV_DIR}/pcgr-linux-64.lock
RUN mamba create -n pcgrr --file ${PCGR_CONDA_ENV_DIR}/pcgrr-linux-64.lock
RUN mamba clean --all --force-pkgs-dirs --yes

FROM quay.io/bioconda/base-glibc-busybox-bash:3.1

COPY --from=0 /opt/mambaforge/envs/ /opt/mambaforge/envs/

ARG PCGR_ENV_NAME="pcgr"
# pcgr env is activated by default
ENV PATH="/opt/mambaforge/envs/${PCGR_ENV_NAME}/bin:${PATH}"
ENV CONDA_PREFIX="/opt/mambaforge/envs/${PCGR_ENV_NAME}"
