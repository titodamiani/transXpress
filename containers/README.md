# SignalP and TargetP

SignalP 6.0 fast and TargetP 2.0 have a DTU licence that forbids sharing them. Each licensee builds an image or a conda env. The main README section "SignalP and TargetP" has the image route.

## Images

`signalp.def` and `targetp.def` build a conda env of fixed packages on `condaforge/miniforge3:26.7.2-0` and copy the DTU package into it. The SignalP image is about 4.6 GiB (the build needs about 12 GB of memory) and the TargetP image about 0.8 GiB. Build commands and the config keys are in the main README. `results/summary/run_info.txt` records the path, size, SIF ID and build time of each image.

## Named conda envs

With `signalp_container: ""` and `targetp_container: ""` the rules use the envs `signalp` and `targetp`.

SignalP 6.0 fast (download from https://services.healthtech.dtu.dk/services/SignalP-6.0/, Downloads):
```bash
conda env create -f workflow/envs/signalp.yaml
conda activate signalp
tar zxvf signalp-6.0h.fast.tar.gz
cd signalp6_fast
pip install signalp-6-package/
SIGNALP_DIR=$(python -c "import signalp; import os; print(os.path.dirname(signalp.__file__))")
cp -r signalp-6-package/models/* $SIGNALP_DIR/model_weights/
conda deactivate
```
Use the python of the conda env.

TargetP 2.0 (download from https://services.healthtech.dtu.dk/services/TargetP-2.0/). The activation script holds the full path of the `targetp-2.0` folder; edit it if you move the folder.
```bash
conda env create -f workflow/envs/targetp.yaml
conda activate targetp
tar zxvf targetp-2.0.Linux.tar.gz
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d" "$CONDA_PREFIX/etc/conda/deactivate.d"
printf 'export _TARGETP_OLD_PATH="$PATH"\nexport PATH="$PATH:%s/targetp-2.0/bin"\n' "$(pwd)" > "$CONDA_PREFIX/etc/conda/activate.d/targetp.sh"
printf 'export PATH="$_TARGETP_OLD_PATH"\nunset _TARGETP_OLD_PATH\n' > "$CONDA_PREFIX/etc/conda/deactivate.d/targetp.sh"
conda deactivate
```
