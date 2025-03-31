from flashtools.compute.data_index import data_index, register_compute_func


@register_compute_func(
    name="r_FLASH", label="$r$", units="cm", cmap="plasma", data_deps=[]
)
def r_FLASH(data, data_yt, **kwargs):
    data["r_FLASH"] = {"data": data_yt["r"].value}
    return data


@register_compute_func(
    name="z_FLASH", label="$Z$", units="cm", cmap="plasma", data_deps=[]
)
def z_FLASH(data, data_yt, **kwargs):
    data["z_FLASH"] = {"data": data_yt["z"].value}
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
    name="depo",
    label="energy deposition",
    units="ergs/g",
    cmap="plasma",
    plot_log10=False,
    data_deps=[],
    data_plot_lims=[-0.05, 0.05],
)
def depo(data, data_yt, **kwargs):
    data["depo"] = {"data": data_yt["depo"].value}
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
    name="tion",
    label=r"$T_i_{\text{FLASH}}$",
    units="C",
    cmap="plasma",
    data_deps=[],
)
def tion(data, data_yt, **kwargs):
    data["tion"] = {"data": data_yt["tion"].value}
    return data


@register_compute_func(
    name="tele",
    label=r"$T_e_{\text{FLASH}}$",
    units="C",
    cmap="plasma",
    data_deps=[],
)
def tele(data, data_yt, **kwargs):
    data["tele"] = {"data": data_yt["tele"].value}
    return data


@register_compute_func(
    name="trad",
    label=r"$T_{\text{rad}_{\text{FLASH}}}$",
    units="C",
    cmap="plasma",
    data_deps=[],
)
def trad(data, data_yt, **kwargs):
    data["trad"] = {"data": data_yt["trad"].value}
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
