"""Compute inputted quantities via recursive function. Inspired by DESC."""

from flashtools.utils import load_2d_data, load_time_series, load_ds, get_FLASH_basis


data_index = {}


def register_compute_func(
    name,
    label,
    units,
    data_deps,
    description=None,
    divergent=False,
    plot_log10=True,
    cmap=None,
    data_plot_lims=None,
    coordinates="rz",
    coordinate_indices=[0, 1],
    basis="rzp",
    **kwargs
):

    deps = {"data": data_deps}

    def _decorator(func):

        d = {
            "fun": func,
            "units": units,
            "label": label,
            "description": description,
            "log": plot_log10,
            "cmap": cmap,
            "data_plot_lims": data_plot_lims,
            "divergent": divergent,
            "coordinates": coordinates,
            "coordinate_indices": coordinate_indices,
            "deps": deps,
            "basis": basis,
            "kwargs": kwargs,
        }
        data_index[name] = d.copy()
        return func

    return _decorator


def compute(
    names,
    # TODO: add object_dir, output_dir
    object_id,
    time_ns=None,
    time_ind=None,
    data_yt=None,
    ts=None,
    ds=None,
    data=None,
    **kwargs
):
    """Have to give some form of an object and some form of a time.
    obj, ts, data_yt, relate to a specific object.
    ds, time_ns, time_ind refer to a specific time within that object
    """
    # could have the dict structure be data[object_id][time_ns][name]
    # or make each new object a class and do obj.compute(name, time_ns)
    if ds is None:
        if ts is None:
            ts = load_time_series(object_id)
        ds = load_ds(ts, time_ns, time_ind)

    # I think this should go above, because you only need ts if you don't have ds and you only need ds
    # if you don't have data_yt
    if data_yt is None:
        data_yt, __ = load_2d_data(object_id, ds=ds)
    if data is None:
        data = {}
    if isinstance(names, str):
        names = [names]
    for name in names:
        if name in data:
            continue
        if len(data_index[name]["deps"]["data"]):
            data = compute(
                data_index[name]["deps"]["data"],
                object_id,
                time_ns=time_ns,
                time_ind=time_ind,
                data_yt=data_yt,
                data=data,
                ts=ts,
                ds=ds,
                **kwargs
            )

        basis = get_FLASH_basis(object_id)
        data = data_index[name]["fun"](data, data_yt, basis=basis, **kwargs)
        data[name].update(
            object_id=object_id, time_ns=ds.current_time.value * 1e9
        )
    return data
