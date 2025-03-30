from .data_index import register_compute_func, data_index
import scipy
import numpy as np
from flashtools.utils import convert_to_eV


@register_compute_func(name="r", label="$r$", units="cm", data_deps=["r_FLASH"])
def r(data, data_yt, **kwargs):
    data["r"] = {"data": np.unique(data["r_FLASH"]["data"][:, :, 0])}
    return data


@register_compute_func(name="z", label="$Z$", units="cm", data_deps=["z_FLASH"])
def z(data, data_yt, **kwargs):
    data["z"] = {"data": np.unique(data["z_FLASH"]["data"][:, :, 0])}
    return data


@register_compute_func(
    name="P_e",
    label=r"$P_\text{e}$",
    units=r"$\text{ergs}/\text{cm}^3$",
    cmap="plasma",
    data_deps=["nele", "T_e"],
    plot_log10=True,
)
def P_e(data, data_yt, **kwargs):
    # convert ev to ergs
    tele = data["T_e"]["data"] * 1.60218e-12
    data["P_e"] = {"data": data["T_e"]["data"] * tele}
    return data


@register_compute_func(
    name="P_i",
    label=r"$P_\text{i}$",
    units=r"$\text{ergs}/\text{cm}^3$",
    cmap="plasma",
    data_deps=["nion", "T_i"],
    plot_log10=True,
)
def P_i(data, data_yt, **kwargs):
    tion = data["T_i"]["data"] * 1.60218e-12
    data["P_i"] = {"data": data["nion"]["data"] * tion}
    return data


@register_compute_func(
    name="P_rad",
    label=r"$P_\text{rad}$",
    units=r"$\text{ergs}/\text{cm}^3$",
    cmap="plasma",
    data_deps=["T_rad"],
    plot_log10=True,
)
def P_rad(data, data_yt, **kwargs):
    a = 4 * scipy.constants.Stefan_Boltzmann / scipy.constants.c
    a *= 1e7  # convert J to ergs
    a *= 1e-6  # convert m^3 to cm^3
    trad = data["T_rad"]["data"] / 8.617e-5  # convert to K
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
    units=data_index["P_tot"]["units"],
    cmap="plasma",
    data_deps=["P_tot"],
    plot_log10=True,
)
def u_tot(data, data_yt, **kwargs):
    data["u_tot"] = {"data": (3 / 2) * data["P_tot"]["data"]}
    return data


# TODO: should have one compute quantity for B and have conversions for
# rzp, xyz in the function
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
    data_plot_lims=[-200, 200]
)
def B_phi(data, data_yt, **kwargs):
    data["B_phi"] = {"data": data["magz"]["data"] * np.sqrt(4 * np.pi) * 1e-3}
    return data


@register_compute_func(
    name="u_mag",
    label=r"$u_{\text{mag}}$",
    units="ergs/cm$^3$",
    description="Internal magnetic energy, see formulary for B units in cgs.",
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
    name="int_u_mag",
    label="$\int r$" + f"{data_index['u_mag']['label']}" + "d$r$",
    units="ergs/cm",
    description="Integrated internal magnetic energy along.",
    cmap="plasma",
    data_deps=["u_mag", "r", "z"],
    divergent=False,
    plot_log10=False,
    coordinates="z"
)
def int_u_mag(data, data_yt, **kwargs):
    z = data["z"]["data"]
    r = data["r"]["data"]
    int_u_mag = np.zeros(len(z))
    for i, z_slice in enumerate(z):
        int_u_mag[i] = np.trapezoid(r * data["u_mag"]["data"][:, i, 0], r)
    data["int_u_mag"] = {"data": int_u_mag}
    return data


@register_compute_func(
    name="T_e",
    label="$T_e$",
    units="eV",
    cmap="plasma",
    data_deps=["tele"],
    data_plot_lims=[0, 4],
)
def T_e(data, data_yt, **kwargs):
    data["T_e"] = {"data": convert_to_eV(data["tele"]["data"])}
    return data


@register_compute_func(
    name="T_i",
    label="$T_i$",
    units="eV",
    cmap="plasma",
    data_deps=["tion"],
)
def T_i(data, data_yt, **kwargs):
    data["T_i"] = {"data": convert_to_eV(data["tion"]["data"])}
    return data


@register_compute_func(
    name="T_rad",
    label=r"$T_{\text{rad}}$",
    units="eV",
    cmap="plasma",
    data_deps=["trad"],
)
def T_rad(data, data_yt, **kwargs):
    data["T_rad"] = {"data": convert_to_eV(data["trad"]["data"])}
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
    r = data["r"]["data"]
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