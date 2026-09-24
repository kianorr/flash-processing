from scipy import constants
import numpy as np
from .data_index import register_compute_func, data_index

@register_compute_func(
    name="m_i",
    label="$m_i$",
    units="kg",
    data_deps=["A"],
)
def ion_mass(data):
    ion_mass = constants.proton_mass * data["A"]
    data["m_i"] = ion_mass
    return data


@register_compute_func(
    name="n_i",
    label="$n_i$",
    units="cm$^{-3}$",
    data_deps=["n_e", "A"],
)
def ion_density(data):
    data["n_i"] = data["n_e"] / data["A"]
    return data


@register_compute_func(
    name="m_e",
    label="$m_e$",
    units="kg",
    data_deps=[],
)
def ion_mass(data):
    data["m_e"] = constants.electron_mass
    return data


@register_compute_func(
    name="q_i",
    label="$q_i$",
    units="C",
    data_deps=["Z"],
)
def ion_charge(data):
    ion_charge = constants.elementary_charge * data["Z"]
    data["q_i"] = ion_charge
    return data


@register_compute_func(
    name="ele_plasma_frequency",
    label="$\omega_{p,e}$",
    units="$1/$s",
    data_deps=["n_e"],
)
def ele_plasma_frequency(data):
    ne_m = data["n_e"] * 100 ** 3
    w = (
        np.sqrt(constants.elementary_charge ** 2 * ne_m / (constants.electron_mass * constants.epsilon_0))
    )
    data["ele_plasma_frequency"] = w
    return data


@register_compute_func(
    name="ion_plasma_frequency",
    label="$\omega_{p,i}$",
    units="$1/$s",
    data_deps=["n_i", "A", "q_i"],
)
def ion_plasma_frequency(data):
    ni_m = data["n_i"] * 100 ** 3
    w = (
        data["q_i"] * np.sqrt(ni_m / (data["A"] * constants.proton_mass * constants.epsilon_0))
    )
    data["ion_plasma_frequency"] = w
    return data


@register_compute_func(
    name="ion_skin_depth",
    label="$\delta_i$",
    units="$\mu$m",
    data_deps=["ion_plasma_frequency"],
)
def ion_skin_depth(data):
    d = constants.c / data["ion_plasma_frequency"]
    data["ion_skin_depth"] = d * 1e6
    return data


@register_compute_func(
    name="electron_skin_depth",
    label="$\delta_e$",
    units="$\mu$m",
    data_deps=["ele_plasma_frequency"],
)
def ele_skin_depth(data):
    d = constants.c / data["ele_plasma_frequency"]
    data["electron_skin_depth"] = d * 1e6
    return data


@register_compute_func(
    name="omega_i",
    label="$\omega_{g,i}$",
    units="s$^{-1}$",
    data_deps=["q_i", "B", "m_i"],
)
def ion_gyrofrequency(data):
    data["omega_i"] = data["q_i"] * data["B"] / data["m_i"]
    return data


@register_compute_func(
    name="omega_e",
    label="$\omega_{g,e}$",
    units="s$^{-1}$",
    data_deps=["B", "m_e"],
)
def ele_gyrofrequency(data):
    data["omega_e"] = constants.elementary_charge * data["B"] / data["m_e"]
    return data


@register_compute_func(
    name="rho_i",
    label="$\\rho_i$",
    units="$\mu$m",
    data_deps=["v_th_i", "omega_i"],
)
def ion_gyroradius(data):
    data["rho_i"] = (data["v_th_i"] / data["omega_i"]) * 1e4
    return data


@register_compute_func(
    name="rho_e",
    label="$\\rho_e$",
    units="$\mu$m",
    data_deps=["v_th_e", "omega_e"],
)
def ele_gyroradius(data):
    data["rho_e"] = (data["v_th_e"] / data["omega_e"]) * 1e4
    return data


@register_compute_func(
    name="B_min",
    label="$B_{\\text{min}}$",
    units="T",
    data_deps=["m_i", "q_i", "v_th_i", "width"],
    description="Minimum magnetic field for trapped ions.",
)
def B_min(data):
    rho_i = data["width"] / 3
    data["B_min"] = data["m_i"] * data["v_th_i"] / (rho_i * data["q_i"])
    return data


@register_compute_func(
    name="B_min_flow",
    label="$B_{\\text{min,flow}}$",
    units="T",
    data_deps=["m_i", "q_i", "v_flow", "width"],
    description="Minimum magnetic field for trapped ions.",
)
def B_min_flow(data):
    rho_i = data["width"] / 3
    data["B_min_flow"] = data["m_i"] * data["v_flow"] / (rho_i * data["q_i"])
    return data


@register_compute_func(
    name="gamma_W",
    label="$\\gamma_{\\text{W}}$",
    units="ns$^{-1}$",
    data_deps=["v_flow", "ion_plasma_frequency"],
    description="Ion weibel growth rate.",
)
def gamma_W(data):
    data["gamma_W"] = (data["v_flow"] * 1e-2) * data["ion_plasma_frequency"] / constants.c * 1e-9
    return data


@register_compute_func(
    name="P_i",
    label="$P_i$",
    units="erg/cm$^3$",
    data_deps=["n_i", "T_i"],
)
def P_i(data):
    eV_to_erg = constants.electron_volt * 1e7
    data["P_i"] = data["n_i"] * data["T_i"] * eV_to_erg
    return data


@register_compute_func(
    name="P_e",
    label="$P_e$",
    units="erg/cm$^3$",
    data_deps=["n_e", "T_e"],
)
def P_e(data):
    eV_to_erg = constants.electron_volt * 1e7
    data["P_e"] = data["n_e"] * data["T_e"] * eV_to_erg
    return data


@register_compute_func(
    name="P_th",
    label="$P_{\\text{th}}$",
    units="erg/cm$^3$",
    data_deps=["P_e", "P_i"],
)
def P_th(data):
    data["P_th"] = data["P_e"] + data["P_i"]
    return data


@register_compute_func(
    name="P_mag",
    label="$P_{\\text{mag}}$",
    units="erg/cm$^3$",
    data_deps=["B"],
)
def P_mag(data):
    data["P_mag"] = data["B"] ** 2 / (2 * constants.mu_0)
    return data


@register_compute_func(
    name="beta",
    label="$\\beta$",
    units="~",
    data_deps=["P_th", "P_mag"],
)
def beta(data):
    data["beta"] = data["P_th"] / data["P_mag"]
    return data


# @register_compute_func(
#     name="coloumb_log",
#     label="$\\Lambda$",
#     units="~",
#     data_deps=[],
# )
# def coloumb_log(data):
#     data["coloumb_log"] = 10
#     return data



@register_compute_func(
    name="kinematic_visc",
    label="$\\nu$",
    units="~",
    data_deps=[],
)
def kinematic_visc(data):
    data["kinematic_visc"] = 2.91e-2
    return data


@register_compute_func(
    name="Re",
    label="$Re$",
    units="~",
    data_deps=["v_flow", "R_length", "kinematic_visc"],
)
def Re(data):
    data["Re"] = data["v_flow"] * 1e-2 * data["R_length"] / data["kinematic_visc"]
    return data


@register_compute_func(
    name="resistivity",
    label="$\\eta$",
    units="~",
    data_deps=[],
)
def resistivity(data):
    data["resistivity"] = 1.14e-7
    return data


@register_compute_func(
    name="Rm",
    label="$Rm$",
    units="~",
    data_deps=["v_flow", "R_length", "resistivity"],
)
def Rm(data):
    # resistivity is in ohm*m
    data["Rm"] = constants.mu_0 * data["v_flow"] * 1e-2 * data["R_length"] / data["resistivity"]
    return data
