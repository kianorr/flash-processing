import scipy
import numpy as np
import os
import scipy.constants
import yt
from flashtools.utils import load_2d_data, load_time_series, convert_to_eV

# TODO: add dim to each variable and make r, z dim=1

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
    coordinates="rzp",
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
    obj,
    time_ns=None,
    time_ind=None,
    data_yt=None,
    ts=None,
    ds=None,
    data=None,
    **kwargs
):
    if ts is None:
        ts = load_time_series(obj)
    if data_yt is None or ds is None:
        data_yt, ds = load_2d_data(obj, time_ns=time_ns, time_index=time_ind, ts=ts)
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
                time_ns=time_ns,
                time_ind=time_ind,
                data_yt=data_yt,
                data=data,
                ts=ts,
                ds=ds,
                **kwargs
            )

        # TODO: don't include these in data_index 
        data_index[name].update(object_id=obj, time_ns=ds.current_time.value * 1e9)
        # data[name] = data_index[name].copy()
        # data[name]["data"] = data_index[name]["fun"](data, data_yt, **kwargs)
        data = data_index[name]["fun"](data, data_yt, **kwargs)
        for key in data_index[name].keys():
            data[name][key] = data_index[name][key]
    return data
