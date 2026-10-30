# The main image of transXpress: the conda env of every rule, built by .github/workflows/image.yml under the tag
# of MAIN_IMAGE in workflow/Snakefile. The env paths come from "snakemake --containerize" on linux-64; each env is
# built from its pin file and post-deploy script. CONTRIBUTING.md, section "The main image".

FROM condaforge/miniforge3:26.7.2-0
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="76eb80a6b4ce21c5e46714879b51208cb1fc9cf3b32dfa362a3c743185f8f135"
LABEL org.opencontainers.image.source="https://github.com/titodamiani/transXpress"

# Conda environment:
#   source: workflow/envs/blast.yaml
#   prefix: /conda-envs/fd730724d7ec7b2ee3441564f1c5d95c
COPY workflow/envs/blast.yaml /conda-envs/fd730724d7ec7b2ee3441564f1c5d95c/environment.yaml
COPY workflow/envs/blast.linux-64.pin.txt /conda-envs/fd730724d7ec7b2ee3441564f1c5d95c/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/fd730724d7ec7b2ee3441564f1c5d95c --file /conda-envs/fd730724d7ec7b2ee3441564f1c5d95c/pin.txt && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/busco.yaml
#   prefix: /conda-envs/c225e36e6ae332ab27564998dd9deb9e
COPY workflow/envs/busco.yaml /conda-envs/c225e36e6ae332ab27564998dd9deb9e/environment.yaml
COPY workflow/envs/busco.linux-64.pin.txt /conda-envs/c225e36e6ae332ab27564998dd9deb9e/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/c225e36e6ae332ab27564998dd9deb9e --file /conda-envs/c225e36e6ae332ab27564998dd9deb9e/pin.txt && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/pfam.yaml
#   prefix: /conda-envs/6132d28153564fd1dbcb90744be549ed
COPY workflow/envs/pfam.yaml /conda-envs/6132d28153564fd1dbcb90744be549ed/environment.yaml
COPY workflow/envs/pfam.linux-64.pin.txt /conda-envs/6132d28153564fd1dbcb90744be549ed/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/6132d28153564fd1dbcb90744be549ed --file /conda-envs/6132d28153564fd1dbcb90744be549ed/pin.txt && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/qc.yaml
#   prefix: /conda-envs/7f6d7640db2982322edcfafa855e6730
COPY workflow/envs/qc.yaml /conda-envs/7f6d7640db2982322edcfafa855e6730/environment.yaml
COPY workflow/envs/qc.linux-64.pin.txt /conda-envs/7f6d7640db2982322edcfafa855e6730/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/7f6d7640db2982322edcfafa855e6730 --file /conda-envs/7f6d7640db2982322edcfafa855e6730/pin.txt && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/rfam.yaml
#   prefix: /conda-envs/001c54639b53066b0221a5cf1540c47c
COPY workflow/envs/rfam.yaml /conda-envs/001c54639b53066b0221a5cf1540c47c/environment.yaml
COPY workflow/envs/rfam.linux-64.pin.txt /conda-envs/001c54639b53066b0221a5cf1540c47c/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/001c54639b53066b0221a5cf1540c47c --file /conda-envs/001c54639b53066b0221a5cf1540c47c/pin.txt && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/rnaspades.yaml
#   prefix: /conda-envs/39f4a4a140ec9685febf193bc07b3778
COPY workflow/envs/rnaspades.yaml /conda-envs/39f4a4a140ec9685febf193bc07b3778/environment.yaml
COPY workflow/envs/rnaspades.linux-64.pin.txt /conda-envs/39f4a4a140ec9685febf193bc07b3778/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/39f4a4a140ec9685febf193bc07b3778 --file /conda-envs/39f4a4a140ec9685febf193bc07b3778/pin.txt && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/tmhmm.yaml
#   prefix: /conda-envs/3774af3c2b2ae9078339f8b44baefa4e
COPY workflow/envs/tmhmm.yaml /conda-envs/3774af3c2b2ae9078339f8b44baefa4e/environment.yaml
COPY workflow/envs/tmhmm.linux-64.pin.txt /conda-envs/3774af3c2b2ae9078339f8b44baefa4e/pin.txt
COPY workflow/envs/tmhmm.post-deploy.sh /conda-envs/3774af3c2b2ae9078339f8b44baefa4e/post-deploy.sh
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/3774af3c2b2ae9078339f8b44baefa4e --file /conda-envs/3774af3c2b2ae9078339f8b44baefa4e/pin.txt && \
    conda run --prefix /conda-envs/3774af3c2b2ae9078339f8b44baefa4e sh /conda-envs/3774af3c2b2ae9078339f8b44baefa4e/post-deploy.sh && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/transdecoder.yaml
#   prefix: /conda-envs/891fd9459bc2c3aebdc0824806fb83b7
COPY workflow/envs/transdecoder.yaml /conda-envs/891fd9459bc2c3aebdc0824806fb83b7/environment.yaml
COPY workflow/envs/transdecoder.linux-64.pin.txt /conda-envs/891fd9459bc2c3aebdc0824806fb83b7/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/891fd9459bc2c3aebdc0824806fb83b7 --file /conda-envs/891fd9459bc2c3aebdc0824806fb83b7/pin.txt && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/trimmomatic.yaml
#   prefix: /conda-envs/64b5bff158d0e06af28fd37216661af9
COPY workflow/envs/trimmomatic.yaml /conda-envs/64b5bff158d0e06af28fd37216661af9/environment.yaml
COPY workflow/envs/trimmomatic.linux-64.pin.txt /conda-envs/64b5bff158d0e06af28fd37216661af9/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/64b5bff158d0e06af28fd37216661af9 --file /conda-envs/64b5bff158d0e06af28fd37216661af9/pin.txt && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/trinity_utils.yaml
#   prefix: /conda-envs/52a57574e85d29e6e4d2485c7d666486
COPY workflow/envs/trinity_utils.yaml /conda-envs/52a57574e85d29e6e4d2485c7d666486/environment.yaml
COPY workflow/envs/trinity_utils.linux-64.pin.txt /conda-envs/52a57574e85d29e6e4d2485c7d666486/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/52a57574e85d29e6e4d2485c7d666486 --file /conda-envs/52a57574e85d29e6e4d2485c7d666486/pin.txt && \
    conda clean --all --yes
