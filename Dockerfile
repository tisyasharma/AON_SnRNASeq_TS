# Container for the AON snRNA-seq pipeline (environment parity, not CI-tested).
# Build:  docker build -t aon-snrnaseq .
# Run:    docker run --rm -it -v "$PWD:/work" aon-snrnaseq snakemake --cores 4 all
FROM mambaorg/micromamba:1.5.10-jammy

WORKDIR /work
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER . /work
CMD ["snakemake", "--cores", "4", "all"]
