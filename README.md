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
```

## compute and plot example
```
from flashtools.plotting import plot_1d, plot_2d
from flashtools.compute import compute, data_index

fig1, ax1 = plt.subplots(1, 1)
fig2, ax2 = plt.subplots(1, 2)
keys = ["nele", "E_dens"]
data = compute(keys, object_dir="/path/to/dir/with/flash_files")
for i, key in enumerate(keys):
    plot_1d(key, data, ax1, slice_of="z", spatial_slice=0.1)
    plot_2d(key, data, ax2[i])
plt.show()
print("the things you can compute are {list(data_index.keys())}")
```
