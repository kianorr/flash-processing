import matplotlib.pyplot as plt
import numpy as np
from flashtools.utils import get_closest, parse_params_file
from flashtools.data_index import compute
from scipy.ndimage import rotate

# TODO: convert units option for input {"units": label, "convert_to": lambda x: x * ...}


def plot_1d(
    name,
    data,
    ax,
    spatial_slice,
    title=False,
    slice_of="z",
    label="",
    return_data=False,
    extent=None,
    conversion=None,
    z_offset=0,
    **kwargs,
):
    """if z, slice should be with the target ending at zero."""
    # TODO: don't need slice_of if you have r, z keywords
    # TODO: make limits actually change data
    obj = data[name]["object_id"]
    target_height = parse_params_file(obj)["sim_targetHeight"]
    log = kwargs.pop("log", data[name]["log"])
    data_1d = data[name]["data"]
    data_1d = data_1d[:, ::-1]
    if conversion is not None:
        data_1d = conversion["convert"](data_1d)
        data_units = conversion["units"]
    else:
        data_units = data[name]['units']

    data_1d = np.log10(data_1d) if log else data_1d

    if "z" not in data:
        data = compute("z", obj, data=data, time_ns=0)
    if "r" not in data:
        data = compute("r", obj, data=data, time_ns=0)
    r = np.unique(data["r"]["data"][:, :, 0])
    z = np.unique(data["z"]["data"][:, :, 0]) + z_offset
    xaxis_name = "r" if slice_of == "z" else "z"

    if slice_of == "z":
        data_1d = data_1d[:, get_closest(z, spatial_slice), 0]
        x_axis = r.copy()
    elif slice_of == "r":
        data_1d = data_1d[get_closest(r, spatial_slice), :, 0]
        x_axis = z.copy()

    # if extent is not None:
    #     for i, lim in enumerate(extent):
    #         if lim is None:
    #             extent[i] = np.min(x_axis) if i % 2 == 0 else np.min(data_1d)
    #             extent[i] = np.max(x_axis) if i % 2 == 1 else np.max(data_1d)
    #             # extent[i] = [r, z][i].min() if i % 2 == 0 else [r, z][i].max()
    #             # extent[i] = [0, 1][i] * [r, z][i].max()
    if extent is not None:
        x_range = (x_axis >= extent[0]) & (x_axis <= extent[1])
    else:
        x_range = slice(None)
        # y_range = (data_1d >= data_1d[x_range]) & (data_1d <= extent[3])
    x_axis = x_axis[x_range]
    data_1d = data_1d[x_range]

    yaxis_label = f"log({data[name]['label']})" if log else data[name]["label"]
    ax.set_xlabel(f"{xaxis_name} [{data[xaxis_name]['units']}]", fontsize=15)
    ax.set_ylabel(f"{yaxis_label} [{data_units}]", fontsize=15)
    if title:
        ax.set_title(
            f"{slice_of} = {spatial_slice} {data[slice_of]['units']}", fontsize=15
        )

    label = f"{label} ({obj})"
    ax.plot(x_axis, data_1d, label=label)

    if return_data:
        return {"x": x_axis, "y": data_1d}


def plot_2d(
    name,
    data,
    ax,
    return_data=False,
    # extent=[-0.1, 0.1, -0.05, 0.45],
    extent=None,
    z_offset=0,
    cbar=False,
    conversion=None,
    **kwargs,
):
    """Plots 2d data.

    Parameters
    ----------
    name: str
        name of data in data dictionary
    data: dict
        dictionary of data
    ax: matplotlib.axes.Axes
        axis to plot on
    extent: list
        [xmin, xmax, ymin, ymax]
    cbar: bool
        whether to include colorbar
    conversion: function
        function to convert data before plotting
    kwargs: dict
        Allowed kwargs are log, data_plot_lims. The rest go to imshow.
    """
    # TODO: generalize to colliding jets

    obj = data[name]["object_id"]
    log = kwargs.pop("log", data[name]["log"])
    # target will always be between (-target_height, 0)
    # target_height = parse_params_file(obj=obj)["sim_targetHeight"]
    if "z" not in data:
        data = compute("z", obj, 0, data=data)
    z = np.unique(data["z"]["data"][:, :, 0])
    if extent is None:
    #     z_range = np.max(z) - np.min(z)
        z += z_offset
        extent=[-0.1, 0.1, np.min(z), np.max(z)]

    data_2d = data[name]["data"][:, :, 0]
    if conversion is not None:
        data_2d = conversion(data_2d)
    data_2d = np.log10(data_2d) if log else data_2d
    # rotate data so that target is at the bottom of the plot
    data_2d = rotate(data_2d, angle=90, reshape=True)
    # data_2d = data_2d[(z >= extent[2]) & (z <= extent[3]), :]
    reflected_data_sign = -1 if data[name]["divergent"] else 1
    data_2d = np.append(reflected_data_sign * np.fliplr(data_2d), data_2d, axis=1)

    data_plot_lims = kwargs.pop(
        "data_plot_lims",
        data[name].get("data_plot_lims", None),
    )
    if data_plot_lims is None:
        data_plot_lims = [np.min(data_2d), np.max(data_2d)]
    p = ax.imshow(
        data_2d,
        extent=extent,
        origin="lower",
        cmap=data[name]["cmap"],
        vmin=data_plot_lims[0],
        vmax=data_plot_lims[1],
        **kwargs,
    )
    ax.grid(False)
    ax.text(
        # -0.09,
        # -0.04,
        extent[0] + 0.01,
        extent[2] + 0.01,
        f"object {obj}",
        ha="left",
        va="bottom",
        color="tab:blue",
        fontsize=8,
    )
    if cbar:
        cbar = plt.gcf().colorbar(p, ax=ax)
        label = f"log({data[name]['label']})" if log else data[name]["label"]
        cbar.set_label(f"{label} [{data[name]['units']}]", fontsize=15)

    if return_data:
        return p, data_2d
    else:
        return p
