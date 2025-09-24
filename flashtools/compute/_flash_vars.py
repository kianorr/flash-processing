from flashtools.compute.data_index import register_compute_func
from flashtools.utils import get_FLASH_basis


@register_compute_func(
    name="first_coord_FLASH", label="$x$ or $r$", units="cm", data_deps=[]
)
def first_coord_FLASH(data, data_yt, **kwargs):
    basis = kwargs.pop("basis", "rzp")
    data["first_coord_FLASH"] = {"data": data_yt[basis[0]].value}
    return data


@register_compute_func(
    name="second_coord_FLASH", label="$y$ or $z$", units="cm", data_deps=[]
)
def second_coord_FLASH(data, data_yt, **kwargs):
    basis = kwargs.pop("basis", "rzp")
    data["second_coord_FLASH"] = {"data": data_yt[basis[1]].value}
    return data

@register_compute_func(
    name="third_coord_FLASH", label="$z$ or None", units="cm", data_deps=[]
)
def second_coord_FLASH(data, data_yt, **kwargs):
    basis = kwargs.pop("basis", "rzp")
    data["third_coord_FLASH"] = {"data": data_yt[basis[2]].value}
    return data

@register_compute_func(name="r_FLASH", label="$r$", units="cm", data_deps=[])
def r_FLASH(data, data_yt, **kwargs):
    data["r_FLASH"] = {"data": data_yt["r"].value}
    return data


@register_compute_func(
    name="z_FLASH", label="$Z$", units="cm", cmap="plasma", data_deps=[]
)
def z_FLASH(data, data_yt, **kwargs):
    data["z_FLASH"] = {"data": data_yt["z"].value}
    return data


@register_compute_func(
    name="x_FLASH", label="$x$", units="cm", cmap="plasma", data_deps=[]
)
def x_FLASH(data, data_yt, **kwargs):
    data["x_FLASH"] = {"data": data_yt["x"].value}
    return data


@register_compute_func(
    name="y_FLASH", label="$y$", units="cm", cmap="plasma", data_deps=[]
)
def y_FLASH(data, data_yt, **kwargs):
    data["y_FLASH"] = {"data": data_yt["y"].value}
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
    name="targ",
    label="target",
    units="~",
    cmap="plasma",
    data_deps=[],
    plot_log10=False,
)
def targ(data, data_yt, **kwargs):
    data["targ"] = {"data": data_yt["targ"].value}
    return data


@register_compute_func(
    name="cham",
    label="chamber",
    units="~",
    cmap="plasma",
    data_deps=[],
    plot_log10=False,
)
def cham(data, data_yt, **kwargs):
    data["cham"] = {"data": data_yt["cham"].value}
    return data


@register_compute_func(
    name="cros",
    label="cross material",
    units="~",
    cmap="plasma",
    data_deps=[],
    plot_log10=False,
)
def cros(data, data_yt, **kwargs):
    data["cros"] = {"data": data_yt["cros"].value}
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
    label=r"$T_{i_{\text{FLASH}}}$",
    units="$\degree$C",
    cmap="plasma",
    data_deps=[],
)
def tion(data, data_yt, **kwargs):
    data["tion"] = {"data": data_yt["tion"].value}
    return data


@register_compute_func(
    name="tele",
    label=r"$T_{e_{\text{FLASH}}}$",
    units="$\degree$C",
    cmap="plasma",
    data_deps=[],
)
def tele(data, data_yt, **kwargs):
    data["tele"] = {"data": data_yt["tele"].value}
    return data


@register_compute_func(
    name="trad",
    label=r"$T_{\text{rad}_{\text{FLASH}}}$",
    units="$\degree$C",
    cmap="plasma",
    data_deps=[],
)
def trad(data, data_yt, **kwargs):
    data["trad"] = {"data": data_yt["trad"].value}
    return data


@register_compute_func(
    name="velx", label="$v_x$", units="cm/s", cmap="plasma", data_deps=[], plot_log10=False
)
def velx(data, data_yt, **kwargs):
    data["velx"] = {"data": data_yt["velx"].value}
    return data


@register_compute_func(
    name="velz", label="$v_z$", units="cm/s", cmap="plasma", data_deps=[], plot_log10=False
)
def velz(data, data_yt, **kwargs):
    data["velz"] = {"data": data_yt["velz"].value}
    return data


@register_compute_func(
    name="vely",
    label="$v_y$",
    units="cm/s",
    cmap="plasma",
    data_deps=[],
    plot_log10=False,
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


@register_compute_func(
    name="eele",
    label="eele",
    units="ergs/g",
    cmap="plasma",
    data_deps=[],
    divergent=False,
    plot_log10=False,
)
def eele(data, data_yt, **kwargs):
    data["eele"] = {"data": data_yt["eele"].value}
    return data


@register_compute_func(
    name="eion",
    label="eion",
    units="ergs/g",
    cmap="plasma",
    data_deps=[],
    divergent=False,
    plot_log10=False,
)
def eion(data, data_yt, **kwargs):
    data["eion"] = {"data": data_yt["eion"].value}
    return data


@register_compute_func(
    name="erad",
    label="erad",
    units="ergs/g",
    cmap="plasma",
    data_deps=[],
    divergent=False,
    plot_log10=False,
)
def erad(data, data_yt, **kwargs):
    data["erad"] = {"data": data_yt["erad"].value}
    return data


def eion(data, data_yt, **kwargs):
    data["eion"] = {"data": data_yt["eion"].value}
    return data


@register_compute_func(
    name="eint",
    label="eint",
    units="ergs/g",
    cmap="plasma",
    data_deps=[],
    divergent=False,
    plot_log10=False,
)
def eint(data, data_yt, **kwargs):
    data["eint"] = {"data": data_yt["eint"].value}
    return data


# TODO: automatically register all FLASH vars that haven't been registered already manually
def register_flash_var(name):
    @register_compute_func(
        name=name,
        label=name,
        units="???",
        cmap="plasma",
        data_deps=[],
        divergent=False,
        plot_log10=False,
    )
    def magp(data, data_yt, **kwargs):
        data["magp"] = {"data": data_yt["magp"].value}
        return data
