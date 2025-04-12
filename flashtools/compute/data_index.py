from flashtools.utils import load_2d_data, load_time_series, load_ds


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
    obj=None,
    time_ns=None,
    object_dir=None,
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
    # could do something where if multiple times are inputted then the dict structure
    # would be data[time][name] and otherwise data[name]
    if ds is None:
        if ts is None:
            ts = load_time_series(obj, object_dir)
        ds = load_ds(ts, time_ns, time_ind)

    # I think this should go above, because you only need ts if you don't have ds and you only need ds
    # if you don't have data_yt
    if data_yt is None:
        data_yt, __ = load_2d_data(obj=obj, object_dir=object_dir, ds=ds)
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
                obj,
                object_dir=object_dir,
                time_ns=time_ns,
                time_ind=time_ind,
                data_yt=data_yt,
                data=data,
                ts=ts,
                ds=ds,
                **kwargs
            )

        # TODO: don't include these in data_index
        data_index[name].update(object_id=obj, object_path=object_dir, time_ns=ds.current_time.value * 1e9)
        data = data_index[name]["fun"](data, data_yt, **kwargs)
        for key in data_index[name].keys():
            data[name][key] = data_index[name][key]
    return data
