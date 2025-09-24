import scipy
import numpy as np
from copy import deepcopy
from .data_index import register_compute_func, data_index
from flashtools.utils import celsius_to_eV


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
    data["P_e"] = {"data": data["nele"]["data"] * tele}
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
    name="e_i",
    label=r"$e_i$",
    units="ergs",
    cmap="plasma",
    data_deps=["P_i"],
    plot_log10=True,
)
def u_i(data, data_yt, **kwargs):
    data["u_i"] = {"data": (3 / 2) * data["P_i"]["data"]}
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


@register_compute_func(
    name="kinetic_energy",
    label=r"$E_{\text{kin}}$",
    units="???",
    cmap="plasma",
    data_deps=["vel_mag"],
    plot_log10=False,
)
def kinetic_energy(data, data_yt, **kwargs):
    data["kinetic_energy"] = {"data": 0.5 * data["vel_mag"]["data"] ** 2}
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
    label=r"$B_z$",
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
    data_plot_lims=[-200, 200],
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
    # magnetic pressure in cgs is B^2 / 8pi
    # where B = B_FLASH * sqrt(4pi)
    FLASH_factor = np.sqrt(4 * np.pi)
    data["u_mag"] = {"data": (data["magz"]["data"] * FLASH_factor) ** 2 / (8 * np.pi)}
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
    data["T_e"] = {"data": celsius_to_eV(data["tele"]["data"])}
    return data


@register_compute_func(
    name="T_i",
    label="$T_i$",
    units="eV",
    cmap="plasma",
    data_deps=["tion"],
    data_plot_lims=[0, 4],
)
def T_i(data, data_yt, **kwargs):
    data["T_i"] = {"data": celsius_to_eV(data["tion"]["data"])}
    return data


@register_compute_func(
    name="T_rad",
    label=r"$T_{\text{rad}}$",
    units="eV",
    cmap="plasma",
    data_deps=["trad"],
)
def T_rad(data, data_yt, **kwargs):
    data["T_rad"] = {"data": celsius_to_eV(data["trad"]["data"])}
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
    data_deps=["first_coord", "second_coord", "velx", "vely"],
    cmap="plasma",
    data_plot_lims=[-1e10, 0.2e10],
    plot_log10=False,
    coordinate_indices=[0, 1],
)
def div_v(data, data_yt, **kwargs):
    basis = kwargs.pop("basis", "rzp")
    first_coord = data["first_coord"]["data"]
    second_coord = data["second_coord"]["data"]
    geometric_factor = first_coord[..., None] if basis == "rzp" else 1

    first_gradient = (
        np.gradient(
            geometric_factor * data["velx"]["data"].squeeze(), first_coord, axis=0
        )
        / geometric_factor
    )
    second_gradient = np.gradient(data["vely"]["data"].squeeze(), second_coord, axis=1)

    div_v = first_gradient + second_gradient

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
    name="vorticity_x",
    label=r"$v$",
    units="cm$^2$/s",
    data_deps=["velz", "vely", "second_coord", "third_coord"],
    cmap="plasma",
    # data_plot_lims=[0, 3500],
    plot_log10=False,
)
def vorticity_x(data, data_yt, **kwargs):
    vz = data["velz"]["data"]
    vy = data["vely"]["data"]
    second_coord = data["second_coord"]["data"]
    third_coord = data["third_coord"]["data"]
    data["vorticity_x"] = {"data": np.gradient(vz, second_coord, axis=1) - np.gradient(vy, third_coord, axis=2)}
    return data


@register_compute_func(
    name="vorticity_y",
    label=r"$v$",
    units="cm$^2$/s",
    data_deps=["velz", "velx", "first_coord", "third_coord"],
    cmap="plasma",
    # data_plot_lims=[0, 3500],
    plot_log10=False,
)
def vorticity_y(data, data_yt, **kwargs):
    vz = data["velz"]["data"]
    vx = data["velx"]["data"]
    first_coord = data["first_coord"]["data"]
    third_coord = data["third_coord"]["data"]
    data["vorticity_y"] = {"data": -(np.gradient(vz, first_coord, axis=0) - np.gradient(vx, third_coord, axis=2))}
    return data


@register_compute_func(
    name="vorticity_z",
    label=r"$v$",
    units="cm$^2$/s",
    data_deps=["vely", "velx", "first_coord", "second_coord"],
    cmap="plasma",
    # data_plot_lims=[0, 3500],
    plot_log10=False,
)
def vorticity_z(data, data_yt, **kwargs):
    vy = data["vely"]["data"]
    vx = data["velx"]["data"]
    first_coord = data["first_coord"]["data"]
    second_coord = data["second_coord"]["data"]
    data["vorticity_z"] = {"data": np.gradient(vy, first_coord, axis=0) - np.gradient(vx, second_coord, axis=1)}
    return data


@register_compute_func(
    name="vorticity_mag",
    label=r"$v$",
    units="cm$^2$/s",
    data_deps=["vorticity_x", "vorticity_y", "vorticity_z"],
    cmap="plasma",
    # data_plot_lims=[0, 3500],
    plot_log10=False,
)
def vorticity_mag(data, data_yt, **kwargs):
    v = np.sqrt(data["vorticity_x"]["data"] ** 2 + data["vorticity_y"]["data"] ** 2 + data["vorticity_z"]["data"] ** 2)
    data["vorticity_mag"] = {"data": v}
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


@register_compute_func(
    name="plasma_beta",
    label="$\\beta$",
    units="~",
    data_deps=["u_mag", "P_tot"],
    cmap="plasma",
    plot_log10=False,
)
def plasma_beta(data, data_yt, **kwargs):
    data["plasma_beta"] = {"data": data["P_tot"]["data"] / data["u_mag"]["data"]}
    return data


@register_compute_func(
    name="beta_percent",
    label=f"{data_index['u_mag']['label']} / ({data_index['u_mag']['label']} + {data_index['P_tot']['label']})",
    units="~",
    data_deps=["u_mag", "P_tot"],
    cmap="plasma",
    plot_log10=False,
)
def beta_percent(data, data_yt, **kwargs):
    total = data["u_mag"]["data"] + data["P_tot"]["data"]
    data["beta_percent"] = {"data": data["u_mag"]["data"] / total}
    return data


@register_compute_func(
    name="j_line",
    label=f"j_line",
    units="~",
    data_deps=["T_e", "nele", "nion"],
    cmap="plasma",
    plot_log10=False,
)
def j_line(data, data_yt, **kwargs):
    j_line = 10 ** -25 * data["nion"]["data"] * data["nele"]["data"] / (data["T_e"]["data"] ** (1 / 2)) * np.exp(-10 / data["T_e"]["data"])
    data["j_line"] = {"data": j_line}
    return data


@register_compute_func(
    name="sound_speed",
    label="$c_s$",
    units="cm/s",
    data_deps=["P_tot", "dens"],
    cmap="plasma",
    plot_log10=False,
)
def sound_speed(data, data_yt, **kwargs):
    c_s = np.sqrt(data["P_tot"]["data"] / data["dens"]["data"])
    data["sound_speed"] = {"data": c_s}
    return data


# TODO: find nice solution to integrating over z instead of r
def integration_1d_helper(first_coord, second_coord, input_data, basis):
    geometric_factor = first_coord if basis == "rzp" else np.ones(len(first_coord))
    # this also works
    # int_1d = np.trapezoid(geometric_factor[..., None] * input_data[:, :, 0], first_coord, axis=0)
    int_1d = np.zeros(len(second_coord))
    for i, z_slice in enumerate(second_coord):
        int_1d[i] = np.trapezoid(geometric_factor * input_data[:, i].squeeze(), first_coord)
    line_element = np.trapezoid(geometric_factor, first_coord)
    int_1d /= line_element
    return int_1d


def integration_2d_helper(first_coord, second_coord, input_data, basis):
    int_1d = integration_1d_helper(first_coord, second_coord, input_data, basis)
    int_2d = np.trapezoid(int_1d, second_coord)
    line_element = (second_coord[-1] - second_coord[0])
    int_2d /= line_element
    return int_2d


@register_compute_func(
    name="int2d_eint",
    label="$\int $" + f"{data_index['eint']['label']}" + "d$A$",
    units=data_index["eint"]["units"],
    description="Volume integrated magnetic internal energy.",
    data_deps=["first_coord", "second_coord", "eint"],
    divergent=False,
    plot_log10=False,
    coordinates="",
    coordinate_indices=[],
)
def int2d_eint(data, data_yt, **kwargs):
    basis = kwargs.pop("basis", "rzp")
    inputs = [data[name]["data"] for name in ["first_coord", "second_coord", "eint"]]
    int_eint = integration_2d_helper(*inputs, basis)
    data["int2d_eint"] = {"data": int_eint}
    return data


@register_compute_func(
    name="int2d_kinetic_energy",
    label="$\int $" + f"{data_index['kinetic_energy']['label']}" + "d$A$",
    units=data_index["kinetic_energy"]["units"],
    description="Volume integrated kinetic energy.",
    data_deps=["first_coord", "second_coord", "kinetic_energy"],
    divergent=False,
    plot_log10=False,
    coordinates="",
    coordinate_indices=[],
)
def int2d_kinetic_energy(data, data_yt, **kwargs):
    basis = kwargs.pop("basis", "rzp")
    inputs = [data[name]["data"] for name in ["first_coord", "second_coord", "kinetic_energy"]]
    int_kinetic_energy = integration_2d_helper(*inputs, basis)
    data["int2d_kinetic_energy"] = {"data": int_kinetic_energy}
    return data


@register_compute_func(
    name="int2d_u_mag",
    label="$\int $" + f"{data_index['u_mag']['label']}" + "d$A$",
    units=data_index["u_mag"]["units"],
    description="Volume integrated magnetic internal energy.",
    data_deps=["first_coord", "second_coord", "u_mag"],
    divergent=False,
    plot_log10=False,
    coordinates="",
    coordinate_indices=[],
)
def int2d_u_mag(data, data_yt, **kwargs):
    basis = kwargs.pop("basis", "rzp")
    inputs = [data[name]["data"] for name in ["first_coord", "second_coord", "u_mag"]]
    int_u_mag = integration_2d_helper(*inputs, basis)
    data["int2d_u_mag"] = {"data": int_u_mag}
    return data


def register_integrations(name):
    @register_compute_func(
        name="int_" + name,
        label=rf"$\left< {(data_index[name]['label']).strip('$')} \right>_\ell$",
        # label="$\int $" + f"{data_index[name]['label']}" + "d$\ell$",
        units=data_index[name]["units"],
        description=f"Integrated {name} over r along z.",
        data_deps=["first_coord", "second_coord", name],
        plot_log10=data_index[name]["log"],
        coordinates="z",
        coordinate_indices=[1]
    )
    def int_(data, data_yt, **kwargs):
        basis = kwargs.pop("basis", "rzp")
        inputs = [data[name]["data"] for name in ["first_coord", "second_coord", name]]
        int_ = integration_1d_helper(*inputs, basis)
        data["int_" + name] = {"data": int_}
        return data
    
di = deepcopy(data_index)
for name in di:
    if len(di[name]["coordinate_indices"]) == 2 and f"int_{name}" not in di:
        register_integrations(name)
