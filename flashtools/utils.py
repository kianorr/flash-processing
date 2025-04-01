import numpy as np
import glob
import re
import os
import yt

yt.set_log_level(50)

# TODO: should probably be a class with the attributes being
# object_dir, object_parent_dir, variables, log_variables


# TODO: global bad bad bad
def set_object_grandparent_dir(new_dir):
    global object_grandparent_dir
    object_grandparent_dir = new_dir
    if not os.path.exists(object_grandparent_dir):
        raise ValueError("Directory does not exist.")


def parse_file(
    obj=None,
    object_dir=None,
    filename="flash.par",
    separator=None,
    break_func=None,
):
    variables = {}
    object_dir = find_path_to_object(object_dir, obj)
    if separator is None:
        separator = ":" if ".log" in filename else "="
    filename = os.path.join(object_dir, filename)
    with open(filename, "r") as file:
        for line in file:
            # Remove comments and strip whitespace
            line = line.split("#")[0].strip()
            if not line:
                continue  # Skip empty lines

            # Match variable assignments (e.g., param=value)
            match = re.match(rf"(.+?)\s*{separator}\s*(.+)", line)
            if match:
                key, value = match.groups()
                # Evaluate value if it's a valid Python literal, else keep as string
                try:
                    value = eval(value)
                except:
                    value = value.strip()
                variables[key] = value

            if break_func is not None and break_func(variables):
                break
    return variables


def parse_log_file(
    obj=None, object_dir=None, filename=None, separator=":", break_func=None
):
    if break_func is None:
        break_func = lambda x: "Comment" in x
    if filename is None:
        filename = find_log_file(obj=obj, object_dir=object_dir)
    variables = parse_file(
        obj=obj,
        object_dir=object_dir,
        filename=filename,
        separator=separator,
        break_func=break_func,
    )
    return variables


def parse_params_file(
    obj=None,
    object_dir=None,
    filename="flash.par",
    separator="=",
    break_func=None,
):
    variables = parse_file(
        obj=obj,
        object_dir=object_dir,
        filename=filename,
        separator=separator,
        break_func=break_func,
    )
    return variables


def get_closest(arr, val):
    """
    Parameters
    ----------
    arr: list or np.ndarray
    val: list, np.ndarray, float

    Returns
    -------
    ind: np.ndarray or number
        if `val` is list, this will be the same shape as `val`
    """
    arr = np.asarray(arr)
    if isinstance(val, list):
        val = np.array(val)
    if not isinstance(val, np.ndarray):
        val = np.array([val])
    ind = np.argmin(np.abs(arr - val[..., None]), axis=1)
    return ind.squeeze()


def find_directory(parent_directory, xxx):
    pattern = f"{parent_directory}/**/object_{xxx}___*"
    matches = glob.glob(pattern)
    directories = [d for d in matches if os.path.isdir(d)]
    if len(directories):
        directories = directories[0]
    else:
        raise ValueError(f"Could not find object {xxx} using {parent_directory}.")

    return directories


def find_path_to_object(object_dir=None, obj=None):
    """
    obj_dir_info: tuple of (parent_directory, object_number)
    object_dir: path to object
    """
    if obj is not None:
        object_dir = find_directory(object_grandparent_dir, obj)
    if object_dir is None:
        raise ValueError("Must provide either an object number or a path to object.")
    return object_dir


def convert_to_eV(temp_C):
    temp_kelvin = temp_C + 273.15
    temp_eV = temp_kelvin * 8.617e-5
    return temp_eV


def find_log_file(obj=None, object_dir=None):
    object_dir = find_path_to_object(object_dir, obj)
    for path in os.listdir(object_dir):
        if ".log" in path:
            with open(os.path.join(object_dir, path), "r") as file:
                if list(file)[0].startswith(" FLASH log file:"):
                    return path


def load_time_series(
    obj=None, object_dir=None, output_dir="output", plot_ext="plt_cnt"
):
    object_dir = find_path_to_object(object_dir, obj)
    variables = parse_params_file(object_dir=object_dir)
    try:
        hdf5_files = os.path.join(
            output_dir, variables["basenm"] + f"hdf5_{plot_ext}_*"
        )
        ts = yt.load(os.path.join(object_dir, hdf5_files))
    except OSError:
        hdf5_files = os.path.join(output_dir, variables["basenm"] + "hdf5_chk_*")
        ts = yt.load(os.path.join(object_dir, hdf5_files))
    return ts


def load_2d_data(
    obj=None,
    object_dir=None,
    ts=None,
    ds=None,
    time_ns=None,
    time_index=None,
    output_dir="output",
):
    object_dir = find_path_to_object(object_dir, obj)
    if ts is None:
        ts = load_time_series(object_dir=object_dir, output_dir=output_dir)
    if time_ns is not None:
        ds = ts.get_by_time((time_ns * 1e-9, "s"))
    elif time_index is not None:
        ds = ts[time_index]
    elif ds is None:
        raise ValueError("Must provide a time.")

    # TODO: these are not time dependent so could be moved to a separate function
    variables = parse_params_file(object_dir=object_dir, filename="flash.par")
    log_variables = parse_log_file(object_dir=object_dir)

    nbx = log_variables["Number x zones"]
    nby = log_variables["Number y zones"]
    data_yt = ds.covering_grid(
        level=variables["lrefine_max"],
        left_edge=np.round(ds.index.grids[0].LeftEdge.value, 2),
        dims=[
            nbx * variables["nblockx"] * 2 ** variables["lrefine_max"],
            nby * variables["nblocky"] * 2 ** variables["lrefine_max"],
            1,
        ],
    )
    # return ds so that it doesn't complain about being temporary
    return data_yt, ds


def compare_dicts(*dicts):
    differences = {}
    all_keys = set()

    for d in dicts:
        all_keys.update(d.keys())

    for key in all_keys:
        values = {i: d.get(key, None) for i, d in enumerate(dicts)}
        unique_values = set(values.values())

        if len(unique_values) > 1:
            differences[key] = values

    return differences


def compare_objects(*objs):
    params = []
    for obj in objs:
        params.append(parse_params_file(obj))

    return compare_dicts(*params)
