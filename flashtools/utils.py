import numpy as np
import glob
import re
import os
import yt

yt.set_log_level(50)


# objects_parent_dir = os.path.abspath("objects")
objects_parent_dir = "/Users/johan/Documents/research/flash_analysis/objects"

def parse_file(
    obj=None,
    path_to_object_dir=None,
    filename="flash.par",
    separator=None,
    break_func=None,
):
    variables = {}
    if obj is not None:
        path_to_object_dir = find_directory(objects_parent_dir, obj)[0]
    if separator is None:
        separator = ":" if ".log" in filename else "="
    filename = os.path.join(path_to_object_dir, filename)
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
    obj=None, path_to_object_dir=None, filename=None, separator=":", break_func=None
):
    if break_func is None:
        break_func = lambda x: "Number x zones" in x and "Number y zones" in x
    if filename is None:
        filename = find_log_file(obj=obj, path_to_object_dir=path_to_object_dir)
    variables = parse_file(
        obj=obj,
        path_to_object_dir=path_to_object_dir,
        filename=filename,
        separator=separator,
        break_func=break_func,
    )
    return variables


def parse_params_file(
    obj=None,
    path_to_object_dir=None,
    filename="flash.par",
    separator="=",
    break_func=None,
):
    variables = parse_file(
        obj=obj,
        path_to_object_dir=path_to_object_dir,
        filename=filename,
        separator=separator,
        break_func=break_func,
    )
    return variables


def get_closest(arr, val):
    ind = np.argmin(np.abs(arr - val))
    return ind


def find_directory(parent_directory, xxx):
    pattern = f"{parent_directory}/**/object_{xxx}___*"
    matches = glob.glob(pattern)
    directories = [d for d in matches if os.path.isdir(d)]

    return directories


def convert_to_eV(temp_C):
    temp_kelvin = temp_C + 273.15
    temp_eV = temp_kelvin * 8.617e-5
    return temp_eV


def find_log_file(obj=None, path_to_object_dir=None):
    # TODO: generalize to par file?
    if obj is not None:
        path_to_object_dir = find_directory(objects_parent_dir, obj)[0]
    for path in os.listdir(path_to_object_dir):
        if ".log" in path:
            with open(os.path.join(path_to_object_dir, path), "r") as file:
                if list(file)[0].startswith(" FLASH log file:"):
                    return path


def load_time_series(obj=None, path_to_object_dir=None, output_dir="output"):
    if obj is not None:
        path_to_object_dir = find_directory(objects_parent_dir, obj)[0]
    variables = parse_params_file(
        obj=obj, path_to_object_dir=path_to_object_dir, filename="flash.par"
    )
    try:
        hdf5_files = os.path.join(output_dir, variables["basenm"] + "hdf5_plt_cnt_*")
    except:
        hdf5_files = os.path.join(output_dir, variables["basenm"] + "hdf5_chk_*")
    ts = yt.load(os.path.join(path_to_object_dir, hdf5_files))
    return ts


def load_2d_data(
    obj=None,
    path_to_object_dir=None,
    ts=None,
    output_dir="output",
    time_ns=None,
    time_index=None,
):
    if ts is None:
        ts = load_time_series(
            obj=obj, path_to_object_dir=path_to_object_dir, output_dir=output_dir
        )
    if obj is not None:
        path_to_object_dir = find_directory(objects_parent_dir, obj)[0]
    if time_ns is not None:
        ds = ts.get_by_time((time_ns * 1e-9, "s"))
    elif time_index is not None:
        ds = ts[time_index]
    else:
        raise ValueError("Must provide a time.")
    
    variables = parse_params_file(
        obj=obj, path_to_object_dir=path_to_object_dir, filename="flash.par"
    )
    log_variables = parse_log_file(obj=obj, path_to_object_dir=path_to_object_dir)
    # hdf5_files = os.path.join(output_dir, variables["basenm"] + "hdf5_chk_*")
    # ts = yt.load(os.path.join(path_to_object_dir, hdf5_files))

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


def compare_dicts(dict1, dict2):
    differences = {}
    all_keys = set(dict1.keys()).union(set(dict2.keys()))

    for key in all_keys:
        val1 = dict1.get(key, None)
        val2 = dict2.get(key, None)

        if val1 != val2:
            differences[key] = {"dict1": val1, "dict2": val2}

    return differences


def compare_objects(obj1, obj2):
    params1 = parse_params_file(obj1)
    params2 = parse_params_file(obj2)

    return compare_dicts(params1, params2)
