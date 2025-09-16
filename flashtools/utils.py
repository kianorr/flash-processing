import numpy as np
import warnings
import glob
import re
import os
import yt

yt.set_log_level(50)


# TODO: global bad bad bad
def set_object_grandparent_dir(new_dir):
    global object_grandparent_dir
    object_grandparent_dir = new_dir
    if not os.path.exists(object_grandparent_dir):
        raise ValueError(f"Directory {new_dir} does not exist.")


# TODO: make separate file for input/output functions
def parse_file(
    object_id,
    filename="flash.par",
    separator=None,
    break_func=None,
):
    variables = {}
    object_dir = find_path_to_object(object_id)
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
    object_id, filename=None, separator=":", break_func=None
):
    if break_func is None:
        break_func = lambda x: "Comment" in x
    if filename is None:
        filename = find_log_file(object_id)
    variables = parse_file(
        object_id=object_id,
        filename=filename,
        separator=separator,
        break_func=break_func,
    )
    return variables


def parse_params_file(
    object_id,
    filename="flash.par",
    separator="=",
    break_func=None,
):
    variables = parse_file(
        object_id,
        filename=filename,
        separator=separator,
        break_func=break_func,
    )
    return variables


def find_directory(parent_directory, xxx):
    pattern = f"{parent_directory}/**/object_{xxx}___*"
    matches = glob.glob(pattern)
    directories = [d for d in matches if os.path.isdir(d)]
    if len(directories):
        directories = directories[0]
    else:
        raise ValueError(f"Could not find object {xxx} using {parent_directory}.")

    return directories


def find_path_to_object(object_id):
    """
    obj_dir_info: tuple of (parent_directory, object_number)
    object_dir: path to object
    """
    if isinstance(object_id, int):
        object_id = find_directory(object_grandparent_dir, object_id)
    elif not isinstance(object_id, str):
        raise ValueError("`object_id` must be a path to the object (str) or number of the object (int).")
    return object_id


def find_log_file(object_id):
    object_dir = find_path_to_object(object_id)
    for path in os.listdir(object_dir):
        if ".log" in path:
            with open(os.path.join(object_dir, path), "r") as file:
                if list(file)[0].startswith(" FLASH log file:"):
                    return path


def load_time_series(
    object_id, output_dir="output", hdf5_extension="plt_cnt"
):
    object_dir = find_path_to_object(object_id)
    variables = parse_params_file(object_id)
    hdf5_files = os.path.join(output_dir, variables["basenm"] + f"hdf5_{hdf5_extension}_*")
    ts = yt.load(os.path.join(object_dir, hdf5_files))
    return ts


def get_times_ns(object_id=None, ts=None):
    if ts is None:
        ts = load_time_series(object_id)
    return [ds.current_time.value * 1e9 for ds in ts]


def load_ds(ts, time_ns=None, time_index=None):
    if time_ns is not None:
        ds = ts.get_by_time((time_ns * 1e-9, "s"))
    elif time_index is not None:
        ds = ts[time_index]
    else:
        raise ValueError("Must provide a time.")

    return ds


def load_2d_data(
    object_id,
    ts=None,
    ds=None,
    time_ns=None,
    time_index=None,
    output_dir="output",
    variables=None,
    log_variables=None,
):
    object_dir = find_path_to_object(object_id)

    if ds is None:
        if ts is None:
            ts = load_time_series(object_id=object_dir, output_dir=output_dir)
        ds = load_ds(ts, time_ns, time_index)

    # TODO: these are not time dependent so could be moved to a separate function
    if variables is None:
        variables = parse_params_file(object_id=object_dir)
    if log_variables is None:
        log_variables = parse_log_file(object_id=object_dir)

    nbx = log_variables["Number x zones"]
    nby = log_variables["Number y zones"]
    nbz = log_variables["Number z zones"]
    z_dim = nbz * variables["nblockz"] * 2 ** variables["lrefine_max"] if "nblockz" in variables else 1
    data_yt = ds.covering_grid(
        level=variables["lrefine_max"],
        left_edge=np.round(ds.index.grids[0].LeftEdge.value, 2),
        dims=[
            nbx * variables["nblockx"] * 2 ** variables["lrefine_max"],
            nby * variables["nblocky"] * 2 ** variables["lrefine_max"],
            z_dim,
        ],
    )
    # return ds so that it doesn't complain about being a weak reference
    return data_yt, ds


# TODO: put these in a separate utils file
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


def compute_resolution(obj, print_res=False):
    variables = parse_params_file(obj)
    log_variables = parse_log_file(obj)
    N_x = (
        log_variables["Number x zones"]
        * variables["nblockx"]
        * 2 ** variables["lrefine_max"]
    )
    N_y = (
        log_variables["Number y zones"]
        * variables["nblocky"]
        * 2 ** variables["lrefine_max"]
    )
    N = N_x * N_y
    return N_x, N_y, N


def convert_to_eV(temp_C):
    warnings.warn("`convert_to_eV` is deprecated. Use `celsius_to_eV` instead.", DeprecationWarning)
    temp_kelvin = temp_C + 273.15
    # multiply by boltzmann constant
    temp_eV = temp_kelvin * 8.617e-5
    return temp_eV


def kelvin_to_eV(temp_kelvin):
    temp_eV = temp_kelvin * 8.617e-5
    return temp_eV


def celsius_to_eV(temp_C):
    return kelvin_to_eV(temp_C + 273.15)


def eV_to_kelvin(temp_eV):
    temp_kelvin = temp_eV / 8.617e-5
    return temp_kelvin


def compare_dicts(*dicts, dict_names=None, exclude=[]):
    if isinstance(exclude, str):
        exclude = [exclude]
    differences = {}
    all_keys = set()

    if dict_names is None:
        dict_names = [i for i in range(len(dicts))]

    for d in dicts:
        all_keys.update(d.keys())

    for key in all_keys:
        if key in exclude:
            continue
        values = {dict_name: d.get(key, None) for dict_name, d in zip(dict_names, dicts)}
        unique_values = set(values.values())

        if len(unique_values) > 1:
            differences[key] = values

    return differences


def compare_objects(*objs, dict_names=None, exclude=[]):
    params = []
    if dict_names is None:
        dict_names = objs
    for obj in objs:
        params.append(parse_params_file(obj))

    return compare_dicts(*params, dict_names=dict_names, exclude=exclude)


def query_files_with_params(filenames, query_params, match_all=True):
    """Return files that contain all or any of the given variable–value pairs."""
    matched_files = []

    for fname in filenames:
        try:
            file_params = parse_file(fname.replace("flash.par", ""))
        except Exception as e:
            print(f"Error reading {fname}: {e}")
            continue

        matches = [
            file_params.get(k) == v
            for k, v in query_params.items()
            if k in file_params
        ]

        if match_all and np.all(matches) and len(matches) == len(query_params):
            matched_files.append(fname)
        elif not match_all and any(matches):
            matched_files.append(fname)

    return matched_files


def query_flashpar_files(query_params, match_all=True):
    """Search for flash.par files that match the given variable–value pairs."""

    # Find all flash.par files in the grandparent directory
    pattern = os.path.join(object_grandparent_dir, "**", "**", "flash.par")
    filenames = glob.glob(pattern)

    # Filter files based on query parameters
    matched_files = query_files_with_params(filenames, query_params, match_all)

    return matched_files


def get_FLASH_basis(obj):
    variables = parse_params_file(obj)
    if variables["geometry"] == "cartesian":
        basis = "xyz"
    elif variables["geometry"] == "cylindrical":
        basis = "rzp"
    return basis
