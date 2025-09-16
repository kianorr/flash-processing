import numpy as np
from .data_index import register_compute_func
from flashtools.utils import get_FLASH_basis


@register_compute_func(
    name="first_coord",
    label="$x$ or $r$",
    units="cm",
    data_deps=["first_coord_FLASH"],
    coordinate_indices=[0],
)
def first_coord(data, data_yt, **kwargs):
    # basis = get_FLASH_basis(data["first_coordinate"]["object_id"])
    # basis = kwargs.pop("basis", "rzp")

    data["first_coord"] = {
        "data": np.unique(data["first_coord_FLASH"]["data"].squeeze())
    }
    return data


@register_compute_func(
    name="second_coord",
    label="$y$ or $z$",
    units="cm",
    data_deps=["second_coord_FLASH"],
    coordinate_indices=[1],
)
def second_coord(data, data_yt, **kwargs):

    data["second_coord"] = {
        "data": np.unique(data["second_coord_FLASH"]["data"].squeeze())
    }
    return data


@register_compute_func(
    name="third_coord",
    label="$z$ or None",
    units="cm",
    data_deps=["third_coord_FLASH"],
    coordinate_indices=[2],
)
def third_coord(data, data_yt, **kwargs):

    data["third_coord"] = {
        "data": np.unique(data["third_coord_FLASH"]["data"].squeeze())
    }
    return data


@register_compute_func(
    name="r", label="$r$", units="cm", data_deps=["r_FLASH"], coordinates="r"
)
def r(data, data_yt, **kwargs):
    data["r"] = {"data": np.unique(data["r_FLASH"]["data"].squeeze())}
    return data


# TODO: bruhh this is gonna be confusing with 3d cartesian
@register_compute_func(
    name="z", label="$Z$", units="cm", data_deps=["z_FLASH"], coordinates="z"
)
def z(data, data_yt, **kwargs):
    data["z"] = {"data": np.unique(data["z_FLASH"]["data"].squeeze())}
    return data


@register_compute_func(
    name="x",
    label="$x$",
    units="cm",
    data_deps=["x_FLASH"],
    coordinates="x",
    basis="xyz",
)
def x(data, data_yt, **kwargs):
    data["x"] = {"data": np.unique(data["x_FLASH"]["data"])}
    return data


@register_compute_func(
    name="y",
    label="$y$",
    units="cm",
    data_deps=["y_FLASH"],
    coordinates="y",
    basis="xyz",
)
def y(data, data_yt, **kwargs):
    data["y"] = {"data": np.unique(data["y_FLASH"]["data"])}
    return data
