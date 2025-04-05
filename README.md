# flashtools

## installation
1. Create an environment

```
conda create -n flash-env python=3.11
conda activate flash-env
```

2. Install the processing code

```
git clone https://github.com/kianorr/flash-processing.git
cd flash-processing
pip install --editable .
pip install -r requirements.txt
```

## compute and plot example
Copy, paste, and replace `object_dir` with the path to a folder containing data, flash.par, and log file.

```
from flashtools.plotting import plot_1d, plot_2d
from flashtools.compute import compute, data_index
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=(8, 4), gridspec_kw={"width_ratios": [0.5, 1]}, dpi=200)

print(f"the things you can compute are {list(data_index.keys())}")
keys = ["nele", "E_dens"]

# dictionary of data
data = compute(keys, time_ns=2.5, object_dir="path/to/object/folder")
print(f"electron density has dimensions (R, Z, Phi) = {np.shape(data['nele']['data'])}")

plot_2d("nele", data, ax[0], cbar=True)
plot_1d("E_dens", data, ax[1], slice_of="z", spatial_slice=0.1)

fig.tight_layout()
plt.show()

```
