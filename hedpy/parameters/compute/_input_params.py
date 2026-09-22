from .data_index import register_compute_func, data_index


@register_compute_func(
    name="system_size",
    label="$L$",
    units="~",
    data_deps=[],
)
def system_size(data):
    return data


@register_compute_func(
    name="A",
    label="$A$",
    units="~",
    data_deps=[],
    description="Atomic number",
)
def A(data):
    return data


@register_compute_func(
    name="A_beta",
    label="$A_{\\beta}$",
    units="~",
    data_deps=[],
    description="Atomic number of second species",
)
def A_beta(data):
    return data


@register_compute_func(
    name="Z",
    label="$Z$",
    units="~",
    data_deps=[],
    description="Charge state.",
)
def Z(data):
    return data


@register_compute_func(
    name="Z_beta",
    label="$Z_{\\beta}$",
    units="~",
    data_deps=[],
    description="Charge state of second species.",
)
def Z_beta(data):
    return data


@register_compute_func(
    name="m_i",
    label="$m_i$",
    units="~",
    data_deps=[],
    description="Ion mass.",
)
def Z_beta(data):
    return data


@register_compute_func(
    name="T_e",
    label="$T_e$",
    units="eV",
    data_deps=[],
    description="Electron temperature.",
)
def T_e(data):
    return data


@register_compute_func(
    name="T_i",
    label="$T_i$",
    units="eV",
    data_deps=[],
    description="Ion temperature.",
)
def T_i(data):
    return data


@register_compute_func(
    name="n_e",
    label="$n_e$",
    units="cm$^{-3}$",
    data_deps=[],
    description="Electron number density.",
)
def n_e(data):
    return data


@register_compute_func(
    name="v_flow",
    label="$v_f$",
    units="cm/s",
    data_deps=[],
    description="Flow velocity from data.",
)
def v_flow(data):
    return data


@register_compute_func(
    name="B",
    label="$B$",
    units="T",
    data_deps=[],
    description="Magnetic field from data.",
)
def B(data):
    return data


@register_compute_func(
    name="dT_e_biermann",
    label="$\\Delta T_e$",
    units="eV",
    data_deps=[],
    description="Change in electron temperature.",
)
def dT_e_biermann(data):
    return data


@register_compute_func(
    name="dne_biermann",
    label="$\\Delta n_e$",
    units="cm$^{-3}$",
    data_deps=[],
    description="Change in electron number density.",
)
def dne_biermann(data):
    return data


@register_compute_func(
    name="dL_ne_biermann",
    label="$\\Delta L_{\\text{ne, Bier}}$",
    units="cm",
    data_deps=[],
    description="Length scale for the Biermann battery effect.",
)
def dL_ne_biermann(data):
    return data


@register_compute_func(
    name="dL_Te_biermann",
    label="$\\Delta L_{\\text{Te, Bier}}$",
    units="cm",
    data_deps=[],
    description="Length scale for the Biermann battery effect.",
)
def dL_Te_biermann(data):
    return data


@register_compute_func(
    name="dt_biermann",
    label="$\\Delta t_{\\text{Biermann}}$",
    units="s",
    data_deps=[],
    description="Time scale for the Biermann battery effect.",
)
def dt_biermann(data):
    return data


@register_compute_func(
    name="R_length",
    label="$l_{Re}$",
    units="cm",
    data_deps=[],
    description="Length scale for Rm and Re.",
)
def l_Re(data):
    return data


@register_compute_func(
    name="coloumb_log",
    label="$\\Lambda$",
    units="~",
    data_deps=[],
    description="Coulomb logarithm.",
)
def coloumb_log(data):
    return data


@register_compute_func(
    name="width",
    label="$\\Delta$",
    units="$\mu$m",
    data_deps=[],
    description="Width of the interaction region.",
)
def width(data):
    return data