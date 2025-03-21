import scipy
import numpy as np
import os
import scipy.constants
import yt
from utils import load_2d_data, load_time_series, convert_to_eV

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
            "deps": deps,
            "kwargs": kwargs,
        }
        data_index[name] = d.copy()
        return func

    return _decorator


@register_compute_func(name="r", label="$r$", units="cm", cmap="plasma", data_deps=[])
def r(data, data_yt, **kwargs):
    data["r"] = {"data": data_yt["r"].value}
    return data


@register_compute_func(name="z", label="$Z$", units="cm", cmap="plasma", data_deps=[])
def z(data, data_yt, **kwargs):
    data["z"] = {"data": data_yt["z"].value}
    return data


@register_compute_func(name="sumy", label="~", units="~", cmap="plasma", data_deps=[])
def sumy(data, data_yt, **kwargs):
    data["sumy"] = {"data": data_yt["sumy"].value}
    return data


@register_compute_func(name="ye", label="ye", units="~", cmap="plasma", data_deps=[])
def ye(data, data_yt, **kwargs):
    data["ye"] = {"data": data_yt["ye"].value}
    return data


@register_compute_func(
    name="shok",
    label="shocks",
    units="~",
    cmap="plasma",
    plot_log10=False,
    data_deps=[],
)
def shok(data, data_yt, **kwargs):
    data["shok"] = {"data": data_yt["shok"].value}
    return data


@register_compute_func(
    name="dens",
    label=r"$\rho$",
    units="g/cm$^3$",
    cmap="plasma",
    data_deps=[],
)
def dens(data, data_yt, **kwargs):
    data["dens"] = {"data": data_yt["dens"].value}
    return data


@register_compute_func(
    name="pres",
    label=r"$P_{\text{FLASH}}$",
    units=r"$\text{dyne}/\text{cm}^2$",
    cmap="plasma",
    data_deps=[],
    plot_log10=False,
)
def pres(data, data_yt, **kwargs):
    data["pres"] = {"data": data_yt["pres"].value}
    return data


@register_compute_func(
    name="P_e",
    label=r"$P_\text{e}$",
    # units=r"$\text{eV}/\text{cm}^3$",
    units=r"$\text{ergs}/\text{cm}^3$",
    cmap="plasma",
    data_deps=["nele", "tele"],
    plot_log10=True,
)
def P_e(data, data_yt, **kwargs):
    # convert ev to ergs
    tele = data["tele"]["data"] * 1.60218e-12
    data["P_e"] = {"data": data["nele"]["data"] * tele}
    return data


@register_compute_func(
    name="P_i",
    label=r"$P_\text{i}$",
    # units=r"$\text{eV}/\text{cm}^3$",
    units=r"$\text{ergs}/\text{cm}^3$",
    cmap="plasma",
    data_deps=["nion", "tion"],
    plot_log10=True,
)
def P_i(data, data_yt, **kwargs):
    tion = data["tion"]["data"] * 1.60218e-12
    data["P_i"] = {"data": data["nion"]["data"] * tion}
    return data


@register_compute_func(
    name="P_rad",
    label=r"$P_\text{rad}$",
    units=r"$\text{ergs}/\text{cm}^3$",
    cmap="plasma",
    data_deps=["trad"],
    plot_log10=True,
)
def P_rad(data, data_yt, **kwargs):
    a = 4 * scipy.constants.Stefan_Boltzmann / scipy.constants.c
    # a *= 6.241509e18 # convert J to eV
    a *= 1e7  # convert J to ergs
    a *= 1e-6  # convert m^3 to cm^3
    trad = data["trad"]["data"] / 8.617e-5  # convert to K
    data["P_rad"] = {"data": a * trad**4 / 3}
    return data


@register_compute_func(
    name="P_tot",
    label=r"$P_{\text{tot}}$",
    units=r"$\text{ergs}/\text{cm}^3$",
    cmap="plasma",
    data_deps=["P_i", "P_e", "P_rad"],
    plot_log10=True,
)
def P_tot(data, data_yt, **kwargs):
    data["P_tot"] = {
        "data": data["P_e"]["data"] + data["P_i"]["data"] + data["P_rad"]["data"]
    }
    return data


@register_compute_func(
    name="u_i",
    label=r"$u_i$",
    # units=r"$\text{eV}/\text{cm}^3$",
    units=data_index["P_i"]["units"],
    cmap="plasma",
    data_deps=["P_i"],
    plot_log10=True,
)
def u_i(data, data_yt, **kwargs):
    data["u_i"] = {"data": (3 / 2) * data["P_i"]["data"]}
    return data


@register_compute_func(
    name="u_e",
    label=r"$u_e$",
    # units=r"$\text{eV}/\text{cm}^3$",
    units=data_index["P_e"]["units"],
    cmap="plasma",
    data_deps=["P_e"],
    plot_log10=True,
)
def u_e(data, data_yt, **kwargs):
    data["u_e"] = {"data": (3 / 2) * data["P_e"]["data"]}
    return data


@register_compute_func(
    name="u_rad",
    label=r"$u_{\text{rad}}$",
    # units=r"$\text{eV}/\text{cm}^3$",
    units=data_index["P_rad"]["units"],
    cmap="plasma",
    data_deps=["P_rad"],
    plot_log10=True,
)
def u_rad(data, data_yt, **kwargs):
    data["u_rad"] = {"data": (3 / 2) * data["P_rad"]["data"]}
    return data


@register_compute_func(
    name="u_tot",
    label=r"$u_{\text{tot}}$",
    # units=r"$\text{eV}/\text{cm}^3$",
    units=data_index["P_tot"]["units"],
    cmap="plasma",
    data_deps=["P_tot"],
    plot_log10=True,
)
def u_tot(data, data_yt, **kwargs):
    data["u_tot"] = {"data": (3 / 2) * data["P_tot"]["data"]}
    return data


@register_compute_func(
    name="velx",
    label="$v_x$",
    units="cm/s",
    cmap="plasma",
    data_deps=[],
)
def velx(data, data_yt, **kwargs):
    data["velx"] = {"data": data_yt["velx"].value}
    return data


@register_compute_func(
    name="vely",
    label="$v_y$",
    units="cm/s",
    cmap="plasma",
    data_deps=[],
)
def vely(data, data_yt, **kwargs):
    data["vely"] = {"data": data_yt["vely"].value}
    return data


@register_compute_func(
    name="magx",
    label="magx",
    description="FLASH output for Br",
    units="G / $\sqrt{4 \pi}$",
    cmap="RdBu",
    data_deps=[],
    divergent=True,
    plot_log10=False,
)
def magx(data, data_yt, **kwargs):
    data["magx"] = {"data": data_yt["magx"].value}
    return data


@register_compute_func(
    name="magy",
    label="magy",
    description="FLASH output for Bz",
    units="G / $\sqrt{4 \pi}$",
    cmap="RdBu",
    data_deps=[],
    divergent=True,
    plot_log10=False,
)
def magy(data, data_yt, **kwargs):
    data["magy"] = {"data": data_yt["magy"].value}
    return data


@register_compute_func(
    name="magz",
    label="magz",
    description="FLASH output for Bphi",
    units="G / $\sqrt{4 \pi}$",
    cmap="RdBu",
    data_deps=[],
    divergent=True,
    plot_log10=False,
)
def magz(data, data_yt, **kwargs):
    data["magz"] = {"data": data_yt["magz"].value}
    return data


@register_compute_func(
    name="B_r",
    label=r"$B_{r}$",
    description="Corrected Br",
    units="kG",
    cmap="RdBu",
    data_deps=["magx"],
    divergent=True,
    plot_log10=False,
)
def B_r(data, data_yt, **kwargs):
    data["B_r"] = {"data": data["magx"]["data"] * np.sqrt(4 * np.pi) * 1e-3}
    return data


@register_compute_func(
    name="B_z",
    label=r"$B_{\phi}$",
    description="Corrected Bz",
    units="kG",
    cmap="RdBu",
    data_deps=["magy"],
    divergent=True,
    plot_log10=False,
)
def B_z(data, data_yt, **kwargs):
    data["B_z"] = {"data": data["magy"]["data"] * np.sqrt(4 * np.pi) * 1e-3}
    return data


@register_compute_func(
    name="B_phi",
    label=r"$B_{\phi}$",
    description="Corrected Bphi",
    units="kG",
    cmap="RdBu",
    data_deps=["magz"],
    divergent=True,
    plot_log10=False,
)
def B_phi(data, data_yt, **kwargs):
    data["B_phi"] = {"data": data["magz"]["data"] * np.sqrt(4 * np.pi) * 1e-3}
    return data


@register_compute_func(
    name="magp",
    label="magp",
    units="???",
    cmap="RdBu",
    data_deps=[],
    divergent=True,
    plot_log10=False,
)
def magp(data, data_yt, **kwargs):
    data["magp"] = {"data": data_yt["magp"].value}
    return data


@register_compute_func(
    name="u_mag",
    label=r"$u_{\text{mag}}$",
    units="ergs/cm$^3$",
    description="Magnetic energy density, see formulary for B units in cgs.",
    cmap="plasma",
    data_deps=["magz"],
    divergent=False,
    plot_log10=False,
)
def u_mag(data, data_yt, **kwargs):
    FLASH_factor = np.sqrt(4 * np.pi)
    data["u_mag"] = {"data": (data["magz"]["data"] * FLASH_factor) ** 2 / (8 * np.pi)}
    return data


@register_compute_func(
    name="tele",
    label="$T_e$",
    units="eV",
    cmap="plasma",
    data_deps=[],
    data_plot_lims=[0, 4],
)
def tele(data, data_yt, **kwargs):
    data["tele"] = {"data": convert_to_eV(data_yt["tele"].value)}
    return data


@register_compute_func(
    name="tion",
    label="$T_i$",
    units="eV",
    cmap="plasma",
    data_deps=[],
)
def tion(data, data_yt, **kwargs):
    data["tion"] = {"data": convert_to_eV(data_yt["tion"].value)}
    return data


@register_compute_func(
    name="trad",
    label=r"$T_{\text{rad}}$",
    units="eV",
    cmap="plasma",
    data_deps=[],
)
def trad(data, data_yt, **kwargs):
    data["trad"] = {"data": convert_to_eV(data_yt["trad"].value)}
    return data


@register_compute_func(
    name="nele",
    label="$n_e$",
    units="cm$^{-3}$",
    data_deps=["ye", "dens"],
    cmap="plasma",
    data_plot_lims=[16, 23],
)
def nele(data, data_yt, **kwargs):
    data["nele"] = {
        "data": data["ye"]["data"] * scipy.constants.N_A * data["dens"]["data"]
    }
    return data


@register_compute_func(
    name="nion",
    label="$n_i$",
    units="cm$^{-3}$",
    data_deps=["sumy", "dens"],
    cmap="plasma",
    data_plot_lims=[16, 23],
)
def nion(data, data_yt, **kwargs):
    data["nion"] = {
        "data": data["sumy"]["data"] * scipy.constants.N_A * data["dens"]["data"]
    }
    return data


@register_compute_func(
    name="div_v",
    label=r"$\nabla \cdot v$",
    units="$1/$s",
    data_deps=["r", "velx", "vely"],
    cmap="plasma",
    data_plot_lims=[-1e10, 0.2e10],
    plot_log10=False,
)
def div_v(data, data_yt, **kwargs):
    r = np.unique(data["r"]["data"][:, :, 0])
    dx = r[1] - r[0]
    div_v = (1 / r) * np.gradient(r * data["velx"]["data"], dx, axis=0) + np.gradient(
        data["vely"]["data"], dx, axis=1
    )
    data["div_v"] = {"data": div_v}
    return data


@register_compute_func(
    name="vel_mag",
    label=r"$v$",
    units="cm/s",
    data_deps=["velx", "vely"],
    cmap="plasma",
    # data_plot_lims=[0, 3500],
    plot_log10=False,
)
def vel_mag(data, data_yt, **kwargs):
    vx = data["velx"]["data"]
    vy = data["vely"]["data"]
    data["vel_mag"] = {"data": np.sqrt(vx**2 + vy**2)}
    return data


@register_compute_func(
    name="E_dens",
    label="$ϵ$",
    units="ergs/cm$^3$",
    data_deps=["dens", "vel_mag"],
    cmap="plasma",
    plot_log10=False,
)
def E_dens(data, data_yt, **kwargs):
    data["E_dens"] = {"data": data["dens"]["data"] * data["vel_mag"]["data"] ** 2 / 2}
    return data


# cm^2/s^2 * g/cm^3 = g/(cm*s^2) = J/cm^3
# g *cm^2/(cm^3*s^2)


def compute(
    names, obj, time_ns=None, time_ind=None, data_yt=None, ts=None, data=None, **kwargs
):
    # TODO: add ability to pass in time_ns as a list
    # especially since all times are loaded in load_2d_data anyway
    if ts is None:
        ts = load_time_series(obj)
    if data_yt is None:
        data_yt, __ = load_2d_data(obj, time_ns=time_ns, time_index=time_ind, ts=ts)
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
                **kwargs
            )

        # TODO: just turn data_index into data? or keep data_index and data separate?
        data_index[name].update(
            object_id=obj,
            time_ns=(
                time_ns
                if time_ns is not None
                else ts[time_ind].current_time.value * 1e9
            ),
        )
        # data[name] = data_index[name].copy()
        # data[name]["data"] = data_index[name]["fun"](data, data_yt, **kwargs)
        data = data_index[name]["fun"](data, data_yt, **kwargs)
        for key in data_index[name].keys():
            data[name][key] = data_index[name][key]
    return data
