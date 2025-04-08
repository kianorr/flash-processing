import matplotlib.pyplot as plt
import numpy as np
from flashtools.utils import get_closest, parse_params_file
from flashtools.compute import compute, data_index
from scipy.ndimage import rotate
import warnings


def plot_1d(
    name,
    data,
    ax,
    # maybe I should be doing
    # xslice=None,
    # yslice=None,
    # zslice=None,
    # or make slice_of a list
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
    data_nd = data[name]["data"].squeeze()
    coordinates = [data_index[name]["coordinates"][i] for i in range(data_nd.ndim)]
    axis_ind_map = {coord: i for i, coord in enumerate(coordinates)}

    if (spatial_slice is None or slice_of is None) and data[name]["data"].ndim > 1:
        raise ValueError("Must provide a slice in r or z for 2D data.")
    if data_nd.ndim == 1 and (spatial_slice is not None or slice_of is not None):
        spatial_slice = None
        slice_of = None

    allowed_kwargs = ["log"]
    for kwarg in kwargs:
        if kwarg not in allowed_kwargs:
            warnings.warn(f"Inputted kwarg {kwarg} not in known kwargs.")

    obj = data[name]["object_id"]
    
    log = kwargs.pop("log", data_index[name]["log"])
    for coord in coordinates:
        if coord not in data:
            obj_path = data[name]["object_path"]
            data = compute(coord, obj, object_dir=obj_path, data=data, time_ns=0)

    # check this so we can reverse array in good conscious
    z_spacing = np.diff(data["z"]["data"])
    if not np.allclose(z_spacing, z_spacing[0], atol=1e-10):
        raise NotImplementedError("z is not uniformly spaced.")

    # reverse array in z (y in xyz?) because target is at the top in FLASH
    axis_to_reverse = list(data_index[name]["basis"])[1]
    data_nd = np.flip(data_nd, axis=axis_ind_map[axis_to_reverse])

    if conversion is not None:
        data_nd = conversion["convert"](data_nd)
        data_units = conversion["units"]
    else:
        data_units = data_index[name]["units"]
    if log:
        data_nd = np.log10(data_nd)

    # could maybe make this more general by popping keys in axis_ind_map
    offset = lambda coord_name: (0 if coord_name == "r" else z_offset)
    xaxis_name = "r" if slice_of == "z" else "z"
    x_axis = data[xaxis_name]["data"] + offset(xaxis_name)

    if spatial_slice is not None:
        spatial_slice_ind = get_closest(
            data[slice_of]["data"] + offset(slice_of), spatial_slice
        )
        s = [slice(None), slice(None)]
        s[axis_ind_map[slice_of]] = slice(spatial_slice_ind, spatial_slice_ind + 1)
        data_1d = data_nd[*s].squeeze()
    else:
        data_1d = data_nd.copy()

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
    ax.plot(x_axis, data_1d, label=label, **kwargs)

    if return_data:
        return {"x": x_axis, "y": data_1d.squeeze()}


def plot_2d(
    name,
    data,
    ax,
    x_range=None,
    y_range=None,
    return_data=False,
    extent=None,
    z_offset=0,
    cbar=False,
    conversion=None,
    title=False,
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

    data_2d = data[name]["data"].squeeze()
    obj = data[name]["object_id"]
    coordinates = [data_index[name]["coordinates"][i] for i in range(data_2d.ndim)]


    allowed_kwargs = ["data_plot_lims", "log"]
    for kwarg in kwargs:
        if kwarg not in allowed_kwargs:
            warnings.warn(f"Inputted kwarg {kwarg} not in known kwargs.")
    log = kwargs.pop("log", data[name]["log"])
    data_plot_lims = kwargs.pop(
        "data_plot_lims",
        data[name].get("data_plot_lims", None),
    )


    for coord in coordinates:
        if coord not in data:
            obj_path = data[name]["object_path"]
            data = compute(coord, obj, 0, obj_path, data=data)
    z = data["z"]["data"].copy()
    z += z_offset
    r = data["r"]["data"].copy()
    r = np.append(-r[::-1], r)


    if conversion is not None:
        data_2d = conversion(data_2d)
    data_2d = np.log10(data_2d) if log else data_2d
    # TODO: make orientation more flexible
    # rotate data so that target is at the bottom of the plot
    data_2d = data_2d.T
    data_2d = data_2d[::-1]
    # mainly for magnetic fields
    reflected_data_sign = -1 if data[name]["divergent"] else 1
    # reflect data across r = 0 axis since data is axisymmetric
    data_2d = np.append(reflected_data_sign * np.fliplr(data_2d), data_2d, axis=1)

    # zoom in to part of the data
    if x_range is None:
        x_range = [np.min(r), np.max(r)]
    if y_range is None:
        y_range = [np.min(z), np.max(z)]
    x_mask = (r >= x_range[0]) & (r <= x_range[1])
    y_mask = (z >= y_range[0]) & (z <= y_range[1])
    data_2d = data_2d[np.ix_(y_mask, x_mask)]

    extent = [*x_range, *y_range]

    if data_plot_lims is None:
        data_plot_lims = [np.min(data_2d), np.max(data_2d)]

    p = ax.imshow(
        data_2d,
        extent=extent,
        origin="lower",
        cmap=data_index[name]["cmap"],
        vmin=data_plot_lims[0],
        vmax=data_plot_lims[1],
        **kwargs,
    )
    if title:
        ax.set_title(f"$t = {data[name]['time_ns']} ns$")
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
        label = f"log({data_index[name]['label']})" if log else data_index[name]["label"]
        cbar.set_label(f"{label} [{data_index[name]['units']}]", fontsize=15)

    if return_data:
        return p, data_2d
    else:
        return p
