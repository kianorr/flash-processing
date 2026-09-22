

data_index = {}

def register_compute_func(
    name,
    label,
    units,
    data_deps,
    description=None,
):

    deps = {"data": data_deps}

    def _decorator(func):

        d = {
            "fun": func,
            "units": units,
            "label": label,
            "deps": deps,
            "description": description,
        }
        data_index[name] = d.copy()
        return func

    return _decorator


def compute(
    names,
    data=None,
):
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
                data=data,
            )

        data = data_index[name]["fun"](data)
    return data


class ComputeData:
    def __init__(self, input_data):
        self.data = input_data

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        self.data[key] = value

    def __contains__(self, key):
        return key in self.data
    
    def compute(self, names):
        if isinstance(names, str):
            names = [names]
        return compute(names, data=self.data)