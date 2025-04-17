import scipy
import numpy as np
from .data_index import register_compute_func, data_index
from flashtools.utils import convert_to_eV, get_FLASH_basis


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


def integration_1d_helper(first_coord, second_coord, input_data, basis):
    geometric_factor = first_coord if basis == "rzp" else 1
    # this also works
    # int_1d = np.trapezoid(geometric_factor[..., None] * input_data[:, :, 0], first_coord, axis=0)
    int_1d = np.zeros(len(second_coord))
    for i, z_slice in enumerate(second_coord):
        int_1d[i] = np.trapezoid(geometric_factor * input_data[:, i, 0], first_coord)
    return int_1d


def integration_2d_helper(first_coord, second_coord, input_data, basis):
    int_1d = integration_1d_helper(first_coord, second_coord, input_data, basis)
    int_2d = np.trapezoid(int_1d, second_coord)
    return int_2d


@register_compute_func(
    # TODO: fix units for cartesian
    name="int_u_mag",
    label="$\int$" + f"{data_index['u_mag']['label']}" + "$d\ell$",
    units="ergs/cm",
    # units=data_index["u_mag"]["units"] + "$d\ell$"
    description="Integrated internal magnetic energy over r along z.",
    cmap="plasma",
    data_deps=["first_coord", "second_coord", "u_mag"],
    divergent=False,
    plot_log10=False,
    coordinates="z",
    coordinate_indices=[1],
)
def int_u_mag(data, data_yt, **kwargs):
    basis = kwargs.pop("basis", "rzp")
    inputs = [data[name]["data"] for name in ["first_coord", "second_coord", "u_mag"]]
    int_u_mag = integration_1d_helper(*inputs, basis)
    data["int_u_mag"] = {"data": int_u_mag}
    return data


@register_compute_func(
    # TODO: fix units for cartesian
    name="int_u_e",
    label="$\int$" + f"{data_index['u_e']['label']}" + "$d\ell$",
    units="ergs/cm",
    # units=data_index["u_mag"]["units"] + "$d\ell$"
    description="Integrated internal electron energy over r along z.",
    cmap="plasma",
    data_deps=["first_coord", "second_coord", "u_e"],
    divergent=False,
    plot_log10=False,
    coordinates="z",
    coordinate_indices=[1],
)
def int_u_e(data, data_yt, **kwargs):
    basis = kwargs.pop("basis", "rzp")
    inputs = [data[name]["data"] for name in ["first_coord", "second_coord", "u_e"]]
    int_u_e = integration_1d_helper(*inputs, basis)
    data["int_u_e"] = {"data": int_u_e}
    return data


@register_compute_func(
    # TODO: fix units for cartesian
    name="int_u_i",
    label="$\int$" + f"{data_index['u_i']['label']}" + "$d\ell$",
    units="ergs/cm",
    # units=data_index["u_mag"]["units"] + "$d\ell$"
    description="Integrated internal electron energy over r along z.",
    cmap="plasma",
    data_deps=["first_coord", "second_coord", "u_i"],
    divergent=False,
    plot_log10=False,
    coordinates="z",
    coordinate_indices=[1],
)
def int_u_i(data, data_yt, **kwargs):
    basis = kwargs.pop("basis", "rzp")
    inputs = [data[name]["data"] for name in ["first_coord", "second_coord", "u_i"]]
    int_u_i = integration_1d_helper(*inputs, basis)
    data["int_u_i"] = {"data": int_u_i}
    return data


@register_compute_func(
    # TODO: fix units for cartesian
    name="int_T_e",
    label="$\int$" + f"{data_index['T_e']['label']}" + "$d\ell$",
    units="eV $\cdot$ cm$^2$",
    # units=data_index["u_mag"]["units"] + "$d\ell$"
    description="Integrated T_e over r along z.",
    cmap="plasma",
    data_deps=["first_coord", "second_coord", "T_e"],
    divergent=False,
    plot_log10=False,
    coordinates="z",
    coordinate_indices=[1],
)
def int_T_e(data, data_yt, **kwargs):
    basis = kwargs.pop("basis", "rzp")
    inputs = [data[name]["data"] for name in ["first_coord", "second_coord", "T_e"]]
    int_T_e = integration_1d_helper(*inputs, basis)
    data["int_T_e"] = {"data": int_T_e}
    return data


@register_compute_func(
    name="int_T_i",
    label="$\int$" + f"{data_index['T_i']['label']}" + "$d\ell$",
    units="eV $\cdot$ cm$^2$",
    description="Integrated T_i over r along z.",
    cmap="plasma",
    data_deps=["first_coord", "second_coord", "T_i"],
    divergent=False,
    plot_log10=False,
    coordinates="z",
    coordinate_indices=[1],
)
def int_T_i(data, data_yt, **kwargs):
    basis = kwargs.pop("basis", "rzp")
    inputs = [data[name]["data"] for name in ["first_coord", "second_coord", "T_i"]]
    int_T_i = integration_1d_helper(*inputs, basis)
    data["int_T_i"] = {"data": int_T_i}
    return data


@register_compute_func(
    # TODO: fix units for cartesian
    name="int_vel_mag",
    label="$\int$" + f"{data_index['vel_mag']['label']}" + "$d\ell$",
    units="cm$^3$/s",
    # units=data_index["u_mag"]["units"] + "$d\ell$"
    description="Integrated flow velocity over r along z.",
    cmap="plasma",
    data_deps=["first_coord", "second_coord", "vel_mag"],
    divergent=False,
    plot_log10=False,
    coordinates="z",
    coordinate_indices=[1],
)
def int_vel_mag(data, data_yt, **kwargs):
    basis = kwargs.pop("basis", "rzp")
    inputs = [data[name]["data"] for name in ["first_coord", "second_coord", "vel_mag"]]
    int_vel_mag = integration_1d_helper(*inputs, basis)
    data["int_vel_mag"] = {"data": int_vel_mag}
    return data


@register_compute_func(
    name="int2d_u_mag",
    label="$\int $" + f"{data_index['u_mag']['label']}" + "d$A$",
    units="ergs",
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


@register_compute_func(
    name="int_plasma_beta",
    label="$\int$" + f"{data_index['plasma_beta']['label']}" + "$d\ell$",
    units="cm$^2$",
    description="Integrated plasma beta over r along z.",
    cmap="plasma",
    data_deps=["first_coord", "second_coord", "plasma_beta"],
    divergent=False,
    plot_log10=False,
    coordinates="z",
    coordinate_indices=[1],
)
def int_plasma_beta(data, data_yt, **kwargs):
    basis = kwargs.pop("basis", "rzp")
    inputs = [data[name]["data"] for name in ["first_coord", "second_coord", "plasma_beta"]]
    int_plasma_beta = integration_1d_helper(*inputs, basis)
    data["int_plasma_beta"] = {"data": int_plasma_beta}
    return data


@register_compute_func(
    name="int_nele",
    label="$\int $" + f"{data_index['nele']['label']}" + "d$\ell$",
    units="cm$^{-1}$",
    description="Integrated electron density energy over r along z.",
    cmap="plasma",
    data_deps=["first_coord", "second_coord", "nele"],
    divergent=False,
    plot_log10=True,
    coordinates="z",
    coordinate_indices=[1],
)
def int_nele(data, data_yt, **kwargs):
    basis = kwargs.pop("basis", "rzp")
    inputs = [data[name]["data"] for name in ["first_coord", "second_coord", "nele"]]
    int_nele = integration_1d_helper(*inputs, basis)
    data["int_nele"] = {"data": int_nele}
    return data


# TODO: implement this with plasmapy units
def register_integrations(name):
    @register_compute_func(
        name="int_" + name,
        label="$\int r$" + f"{data_index[name]['label']}" + "d$r$",
        units="cm$^-1$",
        description=f"Integrated {name} over r along z.",
        cmap="plasma",
        data_deps=[name, "r", "z"],
        divergent=False,
        plot_log10=False,
        coordinates="z",
    )
    # function name might be tricky
    def int_nele(data, data_yt, **kwargs):
        int_name = integration_1d_helper("int_" + name, data)
        data["int_" + name] = {"data": int_name}
        return data
