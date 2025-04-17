import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np
from scipy.ndimage import rotate
import warnings
from flashtools.utils import get_closest, parse_params_file, get_FLASH_basis
from flashtools.compute import compute, data_index


def plot_1d(
    name,
    data,
    ax,
    # maybe I should be doing
    # xslice=None,
    # yslice=None,
    # zslice=None,
    # or make slice_of a list
    slice_of_value=None,
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
    """Plot slice of 2d data or all of 1d data.

    Parameters
    ----------
    name: str
        name of desired plot quantity.
    data: dict
        dictionary of data obtained from compute()
    ax: matplotlib.axes.Axes
        matplotlib axis object. ax will retain whatever is altered in this function.
    """
    data_nd = data[name]["data"].squeeze()
    obj = data[name]["object_id"]
    time = data[name]["time_ns"]

    basis = list(get_FLASH_basis(obj))

    coordinates = [basis[i] for i in data_index[name]["coordinate_indices"]]
    axis_ind_map = {coord: i for i, coord in enumerate(coordinates)}

    should_slice_data = np.all([slice_of_value, slice_of])
    if (slice_of is None) != (slice_of_value is None):
        raise ValueError(
            "Both 'slice_of' and 'slice_of_value' must be provided together or not at all."
        )
    if not should_slice_data and data_nd.ndim > 1:
        raise ValueError("Must provide a slice in r or z for 2D data.")
    elif should_slice_data and data_nd.ndim == 1:
        raise ValueError(
            f"Data is already 1D so {slice_of_value} and {slice_of} are not necessary."
        )

    allowed_kwargs = ["log"]
    for kwarg in kwargs:
        if kwarg not in allowed_kwargs:
            warnings.warn(f"Inputted kwarg {kwarg} not in known kwargs.")

    log = kwargs.pop("log", data_index[name]["log"])
    for coord in coordinates:
        if coord not in data:
            obj_path = data[name]["object_path"]
            data = compute(coord, obj, object_dir=obj_path, data=data, time_ns=0)

    if conversion is not None:
        assert "convert" in conversion and "units" in conversion
        data_nd = conversion["convert"](data_nd)
        data_units = conversion["units"]
    else:
        data_units = data_index[name]["units"]
    if log:
        data_nd = np.log10(data_nd)

    # reverse array in z (y in xyz) because target is at the top in FLASH
    axis_to_reverse = basis[1]
    if axis_to_reverse in axis_ind_map:
        z_spacing = np.diff(data[axis_to_reverse]["data"])
        if not np.allclose(z_spacing, z_spacing[0], atol=1e-10):
            raise NotImplementedError("z is not uniformly spaced.")
        data_nd = np.flip(data_nd, axis=axis_ind_map[axis_to_reverse])

    offset = lambda coord_name: (0 if coord_name == basis[0] else z_offset)

    if should_slice_data:
        slice_of_index = axis_ind_map.pop(slice_of)

        spatial_slice_ind = get_closest(
            data[slice_of]["data"] + offset(slice_of), slice_of_value
        )
        s = [slice(None), slice(None)]
        s[slice_of_index] = slice(spatial_slice_ind, spatial_slice_ind + 1)
        data_1d = data_nd[*s].squeeze()

    else:
        data_1d = data_nd.copy()

    xaxis_name = axis_ind_map.popitem()[0]
    assert not axis_ind_map

    x_axis = data[xaxis_name]["data"] + offset(xaxis_name)

    if x_range is None:
        x_range = [np.min(x_axis), np.max(x_axis)]
    if y_range is None:
        y_range = [np.min(data_1d), np.max(data_1d)]
    x_mask = (x_axis >= x_range[0]) & (x_axis <= x_range[1])
    y_mask = (data_1d >= y_range[0]) & (data_1d <= y_range[1])
    x_axis = x_axis[x_mask & y_mask]
    data_1d = data_1d[x_mask & y_mask]

    yaxis_label = (
        f"log({data_index[name]['label']})" if log else data_index[name]["label"]
    )
    ax.set_xlabel(f"{xaxis_name} [{data_index[xaxis_name]['units']}]", fontsize=15)
    ax.set_ylabel(f"{yaxis_label} [{data_units}]", fontsize=15)
    if title:
        t1 = f"t = {np.round(time, 3)} ns"
        t2 = (
            f" {slice_of} = {slice_of_value} {data_index[slice_of]['units']}"
            if should_slice_data
            else ""
        )
        t = t1 + t2
        ax.set_title(t, fontsize=15)
    obj_text = f"(obj {obj})"
    label = label + rf" $^{{\text{{{obj_text}}}}}$"
    ax.plot(x_axis, data_1d, label=label, **kwargs)

    if return_data:
        return {"x_axis": x_axis, "y_axis": data_1d.squeeze()}


def plot_2d(
    name,
    data,
    ax,
    xaxis_range=None,
    yaxis_range=None,
    yaxis_offset=0,
    conversion=None,
    return_data=False,
    target_loc="bottom",
    cbar=False,
    title=False,
    reflect_across_vertical=False,
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
    cbar: bool
        whether to include colorbar
    conversion: function
        function to convert data before plotting
    kwargs: dict
        Allowed kwargs are log, data_plot_lims. The rest go to imshow.
    """

    data_2d = data[name]["data"].squeeze()
    obj = data[name]["object_id"]
    basis = list(get_FLASH_basis(obj))
    coordinates = [basis[i] for i in data_index[name]["coordinate_indices"]]

    allowed_kwargs = ["data_plot_lims", "log"]
    for kwarg in kwargs:
        if kwarg not in allowed_kwargs:
            warnings.warn(f"Inputted kwarg {kwarg} not in known kwargs.")
    log = kwargs.pop("log", data_index[name]["log"])
    data_plot_lims = kwargs.pop(
        "data_plot_lims",
        data_index[name].get("data_plot_lims", None),
    )

    for coord in coordinates:
        if coord not in data:
            obj_path = data[name]["object_path"]
            data = compute(coord, obj, 0, obj_path, data=data)
    xaxis_name = coordinates[0]
    yaxis_name = coordinates[1]
    y_axis = data[yaxis_name]["data"].copy() + yaxis_offset
    x_axis = data[xaxis_name]["data"].copy()

    if conversion is not None:
        data_2d = conversion(data_2d)
    data_2d = np.log10(data_2d) if log else data_2d
    data_2d = data_2d.T

    assert target_loc in ["top", "bottom", "left", "right"]
    if target_loc == "bottom":
        data_2d = data_2d[::-1]
    elif target_loc == "left" or target_loc == "right":
        raise NotImplementedError("left/right orientation not implented.")

    # reflect data across r = 0 axis since data is axisymmetric
    if xaxis_name == "r" or reflect_across_vertical:
        x_axis = np.append(-x_axis[::-1], x_axis)
        # mainly for magnetic fields
        reflected_data_sign = -1 if data_index[name]["divergent"] else 1
        data_2d = np.append(reflected_data_sign * np.fliplr(data_2d), data_2d, axis=1)

    # zoom in to part of the data
    if xaxis_range is None:
        xaxis_range = [np.min(x_axis), np.max(x_axis)]
    if yaxis_range is None:
        yaxis_range = [np.min(y_axis), np.max(y_axis)]
    x_mask = (x_axis >= xaxis_range[0]) & (x_axis <= xaxis_range[1])
    y_mask = (y_axis >= yaxis_range[0]) & (y_axis <= yaxis_range[1])

    x_axis = x_axis[x_mask]
    y_axis = y_axis[y_mask]
    data_2d = data_2d[np.ix_(y_mask, x_mask)]

    extent = [*xaxis_range, *yaxis_range]

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
        ax.set_title(f"$t = {np.round(data[name]['time_ns'], 3)}$ ns")
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
        label = (
            f"log({data_index[name]['label']})" if log else data_index[name]["label"]
        )
        cbar.set_label(f"{label} [{data_index[name]['units']}]", fontsize=15)

    if return_data:
        return p, {"x": x_axis, "y": y_axis, "data": data_2d}
    else:
        return p


def plot_amr_grid(ds, ax, refinement_filter, widths=None, target_loc="bottom"):
    """Plot amr grid from flash.

    ds: yt.frontends.flash.data_structures.FLASHDataset
        yt dataset
    ax: matplotlib.axes.Axes
        ax to plot on
    refinement_filter: func
        e.g. `lambda level: level != 4` would only show refinement level 4.
        the levels start at zero.
    target_loc: str
        location of target. Default "bottom" flips the grid. 
        put anything else have it normal.
    """


    ymin = ds.domain_left_edge[1]
    ymax = ds.domain_right_edge[1]

    # TODO: make default widths
    # TODO: add option for plotting nxb, nyb

    for grid in ds.index.grids[:]:
        if refinement_filter(grid.Level):
            continue

        left_edge = grid.LeftEdge
        right_edge = grid.RightEdge
        level = grid.Level

        x0 = left_edge[0]
        x1 = right_edge[0]

        y0_orig = left_edge[1]
        y1_orig = right_edge[1]

        # Flip second coordinate (z or y) of the grid
        if target_loc == "bottom":
            y0 = ymax - (y1_orig - ymin)
            y1 = ymax - (y0_orig - ymin)

        rect = patches.Rectangle(
            (x0, y0),
            x1 - x0,
            y1 - y0,
            linewidth=widths[level],
            edgecolor=f"C{level % 10}",
            facecolor="none",
            alpha=0.8,
        )
        ax.add_patch(rect)
