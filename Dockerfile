# The main image of transXpress: the conda env of every rule. .github/workflows/image.yml builds it and
# pushes it under the tag of MAIN_IMAGE in workflow/Snakefile.
#
# The env paths come from "snakemake --containerize" on linux-64, in a T1 run folder (tests/make_run.sh).
# Snakemake finds an env in the image at /conda-envs/<hash>, and the hash covers the env file, its pin file
# for the platform and its post-deploy script. Two edits to the output of snakemake --containerize: a
# fixed tag of the base image, and each env built from its pin file and its post-deploy script, as
# Snakemake does in a run with conda alone. image.yml checks the env paths against snakemake --containerize.

FROM condaforge/miniforge3:26.7.2-0
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="76eb80a6b4ce21c5e46714879b51208cb1fc9cf3b32dfa362a3c743185f8f135"
LABEL org.opencontainers.image.source="https://github.com/titodamiani/transXpress"

# Conda environment:
#   source: workflow/envs/blast.yaml
#   prefix: /conda-envs/fd730724d7ec7b2ee3441564f1c5d95c
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - blast=2.17.0
#     - wget
COPY workflow/envs/blast.yaml /conda-envs/fd730724d7ec7b2ee3441564f1c5d95c/environment.yaml
COPY workflow/envs/blast.linux-64.pin.txt /conda-envs/fd730724d7ec7b2ee3441564f1c5d95c/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/fd730724d7ec7b2ee3441564f1c5d95c --file /conda-envs/fd730724d7ec7b2ee3441564f1c5d95c/pin.txt && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/busco.yaml
#   prefix: /conda-envs/c225e36e6ae332ab27564998dd9deb9e
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - busco=6.1.0
COPY workflow/envs/busco.yaml /conda-envs/c225e36e6ae332ab27564998dd9deb9e/environment.yaml
COPY workflow/envs/busco.linux-64.pin.txt /conda-envs/c225e36e6ae332ab27564998dd9deb9e/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/c225e36e6ae332ab27564998dd9deb9e --file /conda-envs/c225e36e6ae332ab27564998dd9deb9e/pin.txt && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/pfam.yaml
#   prefix: /conda-envs/6132d28153564fd1dbcb90744be549ed
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - hmmer=3.4
#     - wget
COPY workflow/envs/pfam.yaml /conda-envs/6132d28153564fd1dbcb90744be549ed/environment.yaml
COPY workflow/envs/pfam.linux-64.pin.txt /conda-envs/6132d28153564fd1dbcb90744be549ed/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/6132d28153564fd1dbcb90744be549ed --file /conda-envs/6132d28153564fd1dbcb90744be549ed/pin.txt && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/qc.yaml
#   prefix: /conda-envs/7f6d7640db2982322edcfafa855e6730
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - fastqc=0.12.1
#     - multiqc=1.35
COPY workflow/envs/qc.yaml /conda-envs/7f6d7640db2982322edcfafa855e6730/environment.yaml
COPY workflow/envs/qc.linux-64.pin.txt /conda-envs/7f6d7640db2982322edcfafa855e6730/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/7f6d7640db2982322edcfafa855e6730 --file /conda-envs/7f6d7640db2982322edcfafa855e6730/pin.txt && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/rfam.yaml
#   prefix: /conda-envs/001c54639b53066b0221a5cf1540c47c
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - infernal=1.1.5
#     - wget
COPY workflow/envs/rfam.yaml /conda-envs/001c54639b53066b0221a5cf1540c47c/environment.yaml
COPY workflow/envs/rfam.linux-64.pin.txt /conda-envs/001c54639b53066b0221a5cf1540c47c/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/001c54639b53066b0221a5cf1540c47c --file /conda-envs/001c54639b53066b0221a5cf1540c47c/pin.txt && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/rnaspades.yaml
#   prefix: /conda-envs/39f4a4a140ec9685febf193bc07b3778
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - spades=4.3.0
#     - seqkit=2.14.0
COPY workflow/envs/rnaspades.yaml /conda-envs/39f4a4a140ec9685febf193bc07b3778/environment.yaml
COPY workflow/envs/rnaspades.linux-64.pin.txt /conda-envs/39f4a4a140ec9685febf193bc07b3778/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/39f4a4a140ec9685febf193bc07b3778 --file /conda-envs/39f4a4a140ec9685febf193bc07b3778/pin.txt && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/tmhmm.yaml
#   prefix: /conda-envs/3774af3c2b2ae9078339f8b44baefa4e
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     # The same Python as envs/default.yaml, because the tmhmm_parallel script loads the snakemake package
#     # of that env
#     - python=3.11
#     - numpy<1.24            # tmhmm.py 1.2.1 uses np.int, removed in numpy 1.24
#     - cython<3              # builds tmhmm.py from source
#     - gcc
#     - pip
#     - biopython
#     - pip:
#       - tmhmm.py==1.2.1
COPY workflow/envs/tmhmm.yaml /conda-envs/3774af3c2b2ae9078339f8b44baefa4e/environment.yaml
COPY workflow/envs/tmhmm.linux-64.pin.txt /conda-envs/3774af3c2b2ae9078339f8b44baefa4e/pin.txt
COPY workflow/envs/tmhmm.post-deploy.sh /conda-envs/3774af3c2b2ae9078339f8b44baefa4e/post-deploy.sh
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/3774af3c2b2ae9078339f8b44baefa4e --file /conda-envs/3774af3c2b2ae9078339f8b44baefa4e/pin.txt && \
    conda run --prefix /conda-envs/3774af3c2b2ae9078339f8b44baefa4e sh /conda-envs/3774af3c2b2ae9078339f8b44baefa4e/post-deploy.sh && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/transdecoder.yaml
#   prefix: /conda-envs/891fd9459bc2c3aebdc0824806fb83b7
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     # The rules run the two TransDecoder stages, TransDecoder.LongOrfs and TransDecoder.Predict.
#     # The package keeps them in util/ next to the TransDecoder wrapper, outside bin/, and the
#     # rules find them from the path of the wrapper.
#     # Build 0 of 6.0.0 links no TransDecoder wrapper into bin/, so the pin asks for build 1 or later.
#     - "transdecoder[version='6.0.0',build_number='>=1']"
COPY workflow/envs/transdecoder.yaml /conda-envs/891fd9459bc2c3aebdc0824806fb83b7/environment.yaml
COPY workflow/envs/transdecoder.linux-64.pin.txt /conda-envs/891fd9459bc2c3aebdc0824806fb83b7/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/891fd9459bc2c3aebdc0824806fb83b7 --file /conda-envs/891fd9459bc2c3aebdc0824806fb83b7/pin.txt && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/trimmomatic.yaml
#   prefix: /conda-envs/64b5bff158d0e06af28fd37216661af9
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - trimmomatic=0.41
COPY workflow/envs/trimmomatic.yaml /conda-envs/64b5bff158d0e06af28fd37216661af9/environment.yaml
COPY workflow/envs/trimmomatic.linux-64.pin.txt /conda-envs/64b5bff158d0e06af28fd37216661af9/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/64b5bff158d0e06af28fd37216661af9 --file /conda-envs/64b5bff158d0e06af28fd37216661af9/pin.txt && \
    conda clean --all --yes

# Conda environment:
#   source: workflow/envs/trinity_utils.yaml
#   prefix: /conda-envs/52a57574e85d29e6e4d2485c7d666486
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - trinity=2.13.2=hea94271_3
#     # Trinity calls salmon with --validateMappings, which salmon 2 removed.
#     # The bioconda trinity recipe sets no bound, so salmon must be pinned here.
#     - salmon=1.10.3
#     - bioconductor-edger=3.36.0=r41hc247a5b_2
#     - kallisto=0.48.0=h15996b6_2
#     - r=4.1=r41hd8ed1ab_1006
#     - r-tidyverse=1.3.2=r41hc72bb7e_1
#     - bowtie2=2.5.0=py310h8d7afc0_0
#     - samtools=1.16.1=h6899075_1
#     - rsem=1.3.3=pl5321ha04fe3b_5
COPY workflow/envs/trinity_utils.yaml /conda-envs/52a57574e85d29e6e4d2485c7d666486/environment.yaml
COPY workflow/envs/trinity_utils.linux-64.pin.txt /conda-envs/52a57574e85d29e6e4d2485c7d666486/pin.txt
RUN conda create --quiet --yes --no-default-packages --prefix /conda-envs/52a57574e85d29e6e4d2485c7d666486 --file /conda-envs/52a57574e85d29e6e4d2485c7d666486/pin.txt && \
    conda clean --all --yes
