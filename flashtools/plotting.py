import matplotlib.pyplot as plt
import numpy as np
from flashtools.utils import get_closest, parse_params_file
from flashtools.compute import compute, data_index
from scipy.ndimage import rotate
import warnings

# TODO: convert units option for input {"units": label, "convert_to": lambda x: x * ...}


def plot_1d(
    name,
    data,
    ax,
    # maybe I should be doing
    # xslice=None,
    # yslice=None,
    # zslice=None,
    spatial_slice=None,
    slice_of=None,
    x_range=None,
    y_range=None,
    z_offset=0,
    title=False,
    label="",
    conversion=None,
    return_data=False,
    **kwargs,
):
    """if z, slice should be with the target ending at zero."""
    data_1d = data[name]["data"].squeeze()
    coordinates = [data_index[name]["coordinates"][i] for i in range(data_1d.ndim)]
    # axis_ind_map = {data_index[name]["basis"][i]: i for i in range(data_1d.ndim)}
    axis_ind_map = {"r": 0, "z": 1, "p": 2}

    if (spatial_slice is None or slice_of is None) and data[name]["data"].ndim > 1:
        raise ValueError("Must provide a slice in r or z for 2D data.")

    allowed_kwargs = ["log"]
    for kwarg in kwargs:
        if kwarg not in allowed_kwargs:
            warnings.warn(f"Inputted kwarg {kwarg} not in known kwargs.")
    
    obj = data[name]["object_id"]
    log = kwargs.pop("log", data_index[name]["log"])
    for coord in coordinates:
        if coord not in data:
            data = compute(coord, obj, data=data, time_ns=0)

    # check this so we can reverse array in good conscious
    z_spacing = np.diff(data["z"]["data"])
    if not np.allclose(z_spacing, z_spacing[0], atol=1e-10):
        raise NotImplementedError("z is not uniformly spaced.")
    
    # reverse array in Z because target is at the top in FLASH
    s = np.array([slice(None), slice(None, None, -1), slice(None)])
    s = s[[axis_ind_map[coordinates[i]] for i in range(data_1d.ndim)]]
    data_1d = data_1d[*s]

    if conversion is not None:
        data_1d = conversion["convert"](data_1d)
        data_units = conversion["units"]
    else:
        data_units = data_index[name]['units']
    if log:
        data_1d = np.log10(data_1d)

    # could maybe make this more general by popping keys in axis_ind_map
    offset = lambda coord_name: (0 if coord_name == "r" else z_offset)
    xaxis_name = "r" if slice_of == "z" else "z"
    x_axis = data[xaxis_name]["data"] + offset(xaxis_name)


    if spatial_slice is not None:
        spatial_slice_ind = get_closest(data[slice_of]["data"] + offset(slice_of), spatial_slice)
        s = [slice(None), slice(None)]
        s[axis_ind_map[slice_of]] = slice(spatial_slice_ind, spatial_slice_ind + 1)
        data_1d = data_1d[*s].squeeze()


    if x_range is None:
        x_range = [np.min(x_axis), np.max(x_axis)]
    if y_range is None:
        y_range = [np.min(data_1d), np.max(data_1d)]
    x_mask = (x_axis >= x_range[0]) & (x_axis <= x_range[1])
    y_mask = (data_1d >= y_range[0]) & (data_1d <= y_range[1])
    x_axis = x_axis[x_mask & y_mask]
    data_1d = data_1d[x_mask & y_mask]

    yaxis_label = f"log({data[name]['label']})" if log else data[name]["label"]
    ax.set_xlabel(f"{xaxis_name} [{data[xaxis_name]['units']}]", fontsize=15)
    ax.set_ylabel(f"{yaxis_label} [{data_units}]", fontsize=15)
    if title:
        ax.set_title(
            f"{slice_of} = {spatial_slice} {data_index[slice_of]['units']}", fontsize=15
        )
    obj_text = f"(obj {obj})"
    label = label + rf" $^{{\text{{{obj_text}}}}}$"
    ax.plot(x_axis, data_1d, label=label)

    if return_data:
        return {"x": x_axis, "y": data_1d.squeeze()}


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

    allowed_kwargs = ["data_plot_lims", "log"]
    for kwarg in kwargs:
        if kwarg not in allowed_kwargs:
            warnings.warn(f"Inputted kwarg {kwarg} not in known kwargs.")

    obj = data[name]["object_id"]
    log = kwargs.pop("log", data[name]["log"])
    if "z" not in data:
        data = compute("z", obj, 0, data=data)
    z = data["z"]["data"]
    if extent is None:
        z += z_offset
        extent=[-0.1, 0.1, np.min(z), np.max(z)]

    data_2d = data[name]["data"].squeeze()
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
