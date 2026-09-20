from scipy import constants
import numpy as np
from .data_index import register_compute_func, data_index


@register_compute_func(
    name="v_th_i",
    label="$v_{\\text{th}, i}$",
    units="cm/s",
    data_deps=["T_i", "m_i"],
)
def thermal_velocity(data):
    v = np.sqrt(
        constants.electron_volt * data["T_i"] / data["m_i"]
    )
    data["v_th_i"] = v * 1e2
    return data


@register_compute_func(
    name="v_th_e",
    label="$v_{\\text{th}, e}$",
    units="cm/s",
    data_deps=["T_e"],
)
def ele_thermal_velocity(data):
    v = np.sqrt(
        constants.electron_volt * data["T_e"] / constants.electron_mass
    )
    data["v_th_e"] = v * 1e2
    return data


@register_compute_func(
    name="v_th_ie",
    label="$v_{\\text{th}, ie}$",
    units="cm/s",
    data_deps=["v_th_i", "v_th_e"],
)
def ion_ele_thermal_velocity(data):
    v = np.sqrt(data["v_th_i"] ** 2 + data["v_th_e"] ** 2)
    # v = np.sqrt(
    #     constants.electron_volt * data["T_e"] / constants.electron_mass
    # )
    data["v_th_ie"] = v
    return data


@register_compute_func(
    name="T_i_flow",
    label="$T_{i, \\text{thermalized}}$",
    units="keV",
    data_deps=["v_flow", "m_i"],
    description="T_i corresponding to flow -> thermalization"
)
def T_i_flow(data):
    T_i = (data["v_flow"] * 1e-2) ** 2 * data["m_i"] / (2 * constants.electron_volt)
    data["T_i_flow"] = T_i * 1e-3
    return data


@register_compute_func(
    name="ion_collision_freq",
    label="$\\nu_i$",
    units="1/s",
    data_deps=["n_i", "T_i", "Z", "m_i"],
    description="from formulary, page 28"
)
def ion_collision_freq(data):
    mu = data["m_i"] / constants.proton_mass
    data["ion_collision_freq"] = (
        4.8e-8 * data["Z"] ** 4 * data["n_i"] * data["coloumb_log"] / (np.sqrt(mu) * data["T_i"] ** (3/2))
    )
    return data


@register_compute_func(
    name="ele_collision_freq",
    label="$\\nu_e$",
    units="1/s",
    data_deps=["n_e", "T_e"],
)
def ele_collision_freq(data):
    data["ele_collision_freq"] = (
        2.91e-6 * data["n_e"] * data["coloumb_log"] * data["T_e"] ** (-3/2)
    )
    return data


@register_compute_func(
    name="ion_mfp",
    label="$\\lambda_{\\text{th},i}$",
    units="$\\mu$m",
    data_deps=["ion_collision_freq", "v_th_i"],
)
def ion_mfp(data):
    data["ion_mfp"] = data["v_th_i"] / data["ion_collision_freq"]
    data["ion_mfp"] *= 1e4
    return data


@register_compute_func(
    name="ele_mfp",
    label="$\\lambda_{\\text{th},e}$",
    units="$\\mu$m",
    data_deps=["ele_collision_freq", "v_th_e"],
)
def ele_mfp(data):
    data["ele_mfp"] = data["v_th_e"] / data["ele_collision_freq"]
    # converting cm to microns
    data["ele_mfp"] *= 1e4
    return data


@register_compute_func(
    name="ion-ion_mfp",
    label="$\\lambda_{ii}$",
    units="$\\mu$m",
    data_deps=["A", "A_beta", "Z", "Z_beta", "n_i", "v_flow", "coloumb_log"],
    description="Fiuza 2020 supplementary; equation 2"
)
def ion_ion_mfp(data):
    n_beta = data["n_i"] / 1e19
    v_flow = data["v_flow"] / 1e8
    numerator = 670 * data["A"] ** 2 * data["A_beta"] ** 2 * v_flow ** 4
    denominator = data["Z"] ** 2 * data["Z_beta"] ** 2 * (data["A"] + data["A_beta"]) ** 2 * n_beta * data["coloumb_log"]
    # return numerator / denominator
    data["ion-ion_mfp"] = numerator / denominator * 1e4
    return data


@register_compute_func(
    name="ion-ele_mfp",
    label="$\\lambda_{ie}$",
    units="$\\mu$m",
    data_deps=["Z", "n_i", "T_e", "coloumb_log"],
    description="Fiuza 2020 supplementary; equation 3"
)
def ion_ele_mfp(data):
    n_beta = data["n_i"] / 1e19
    numerator = 0.4 * (data["T_e"] / 1000) ** 2
    denominator = data["Z"] ** 2 * n_beta * data["coloumb_log"]
    data["ion-ele_mfp"] = numerator / denominator * 1e4
    return data


@register_compute_func(
    name="ele_heating_time",
    label="$\\tau_{ie}$",
    units="$ns$",
    data_deps=["Z", "n_e", "T_e", "coloumb_log"],
    description="Equation 5 of Fiuza 2020 supplementary"
)
def ion_ele_mfp(data):
    n_beta = data["n_i"] / 1e19
    v_flow = data["v_flow"] / 1e8
    T_e = data["T_e"] / 1000
    t = 38 * T_e ** (5 / 2) / (data["Z"] * n_beta * v_flow ** 2 * data["coloumb_log"])
    data["ele_heating_time"] = t
    return data