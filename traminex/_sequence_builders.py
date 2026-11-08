from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
import pandas as pd

from ._categorical_order import r_seqstatl_order
from .sequences import STSList

SEQ_KWARGS = frozenset(
    {
        "labels",
        "cpal",
        "xtlab",
        "xtstep",
        "tick_last",
        "nr",
        "void",
        "weights",
        "start",
        "ids",
        "missing_color",
    }
)


def validate_seq_kwargs(kwargs: dict, *, func_name: str) -> None:
    unknown = sorted(set(kwargs) - SEQ_KWARGS)
    if unknown:
        raise TypeError(
            f"Unsupported kwargs for {func_name}: {unknown}. Allowed kwargs: {sorted(SEQ_KWARGS)}"
        )
    nr, void = kwargs.get("nr", "*"), kwargs.get("void", "%")
    if nr is not None and void is not None and str(nr) == str(void):
        # seqdef's refusal: R would read a missing cell as void (seqdef ledger row,
        # TX-AUD-20261010-035); the builders respelled one marker (TX-REP-20261019-012).
        raise ValueError("'nr' and 'void' must be distinct character markers")


def infer_alphabet(X: np.ndarray) -> list[str]:
    """The distinct states of ``X`` in R's ``seqdef`` order (``seqstatl``: numbers by
    value when every state is an R number, R's collation otherwise; owner decision D2)."""
    flat = np.asarray(X, dtype=object).astype(str).ravel(order="C")
    return r_seqstatl_order([str(state) for state in pd.unique(flat).tolist()])


def _seqdef_default_missing(missing: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Split missing cells as TraMineR's ``seqdef`` defaults do.

    ``left = NA`` and ``gaps = NA`` leave leading and internal missing cells as
    nonresponse; ``right = "DEL"`` makes the trailing ones -- every cell of a
    row with no observed state -- void (``SOURCES/TraMineR/R/seqdef.R:7``).
    Returns the nonresponse and void masks.
    """
    trailing = np.zeros(missing.shape, dtype=bool)
    for row in range(missing.shape[0]):
        observed = np.flatnonzero(~missing[row])
        trailing[row, (int(observed[-1]) + 1 if observed.size else 0) :] = True
    return missing & ~trailing, missing & trailing


# seqdef's default palette for the real states (seqdef.R: RColorBrewer's Accent up to 8
# states, Set3 up to 12, else colorspace's qualitative_hcl), shared by seqdef and the builders.
_ACCENT = [
    "#7FC97F",
    "#BEAED4",
    "#FDC086",
    "#FFFF99",
    "#386CB0",
    "#F0027F",
    "#BF5B17",
    "#666666",
]
_SET3 = [
    "#8DD3C7",
    "#FFFFB3",
    "#BEBADA",
    "#FB8072",
    "#80B1D3",
    "#FDB462",
    "#B3DE69",
    "#FCCDE5",
    "#D9D9D9",
    "#BC80BD",
    "#CCEBC5",
    "#FFED6F",
]


def _linear_to_srgb(value: float) -> float:
    if value > 0.00304:
        return 1.055 * value ** (1.0 / 2.4) - 0.055
    return 12.92 * value


def _qualitative_hcl_color(hue: float, *, lightness: float = 80, chroma: float = 50) -> str:
    """Convert one colorspace Set-3 polar-LUV color to an sRGB hex code."""
    angle = math.radians(hue)
    u_star = chroma * math.cos(angle)
    v_star = chroma * math.sin(angle)

    xn, yn, zn = 95.047, 100.0, 108.883
    white_denom = xn + 15.0 * yn + 3.0 * zn
    u_white = 4.0 * xn / white_denom
    v_white = 9.0 * yn / white_denom
    y = yn * ((lightness + 16.0) / 116.0) ** 3
    u_prime = u_star / (13.0 * lightness) + u_white
    v_prime = v_star / (13.0 * lightness) + v_white
    x = y * 9.0 * u_prime / (4.0 * v_prime)
    z = y * (12.0 - 3.0 * u_prime - 20.0 * v_prime) / (4.0 * v_prime)

    x /= 100.0
    y /= 100.0
    z /= 100.0
    linear = (
        3.240479 * x - 1.537150 * y - 0.498535 * z,
        -0.969256 * x + 1.875992 * y + 0.041556 * z,
        0.055648 * x - 0.204043 * y + 1.057311 * z,
    )
    channels = [min(max(_linear_to_srgb(channel), 0.0), 1.0) for channel in linear]
    rgb = [int(math.floor(channel * 255.0 + 0.5)) for channel in channels]
    return "#" + "".join(f"{channel:02X}" for channel in rgb)


def _default_seqdef_palette(size: int) -> list[str]:
    if size <= 0:
        return []
    if size <= 8:
        return _ACCENT[:size]
    if size <= 12:
        return _SET3[:size]
    return [_qualitative_hcl_color(10.0 + 360.0 * index / size) for index in range(size)]


def encode_state_matrix(
    X: np.ndarray,
    *,
    nr: object = "*",
    void: object = "%",
) -> tuple[np.ndarray, list[str], dict[str, int], dict[str, str | None]]:
    """Encode a state matrix; missing cells follow ``seqdef``'s defaults.

    Observed cells keep their text and define the alphabet, in R's ``seqdef``
    order (:func:`infer_alphabet`). A missing cell
    (``None``, NaN, ``pd.NA``) becomes the ``nr`` or ``void`` marker by
    :func:`_seqdef_default_missing`, and a marker that is used is appended after
    the real states, nonresponse first, as ``seqdef`` does (audit S2-05). A real
    state spelled like a marker stays a state: the marker then gets another
    spelling, as on seqdef's SPS and SPELL routes (TX-AUD-20261019-005). Returns
    the codes, the alphabet, the codebook and the marker texts used (``None`` for
    a marker no cell needs).
    """
    from .seqformat import _collision_free_marker, _r_state_text

    X_obj = np.asarray(X, dtype=object)
    missing = np.asarray(pd.isna(X_obj), dtype=bool)
    # Each observed cell as R prints it, as seqdef reads it: 1.0 is "1"
    # (TX-AUD-20261010-073).
    X_str = np.empty(X_obj.shape, dtype=object)
    for index, value in np.ndenumerate(X_obj):
        X_str[index] = "" if missing[index] else _r_state_text(value)
    alphabet = infer_alphabet(X_str[~missing])
    used: dict[str, str | None] = {"nr": None, "void": None}
    if missing.any():
        for cells, marker, name, role in zip(
            _seqdef_default_missing(missing),
            (nr, void),
            ("nr", "void"),
            ("NR", "VOID"),
            strict=True,
        ):
            if not cells.any():
                continue
            if marker is None:
                raise ValueError(
                    f"the state table has missing cells that need the {name} marker, "
                    f"but {name}=None disables it"
                )
            text = _collision_free_marker(str(marker), alphabet, role)
            X_str[cells] = text
            alphabet.append(text)
            used[name] = text
    code = {state: i for i, state in enumerate(alphabet)}
    seq = np.vectorize(code.__getitem__, otypes=[int])(X_str)
    return seq, alphabet, code, used


def build_sts_from_state_matrix(X: np.ndarray, *, kwargs: dict | None = None) -> STSList:
    """An STSList from a state table; ``labels`` and ``cpal`` pair with the states in
    the alphabet's order, which is R's ``seqdef`` order, as in :func:`traminex.seqdef`.
    Marker codes are set for the markers the table uses only, so a real state spelled
    like an unused marker is not read as one."""
    options = {} if kwargs is None else kwargs
    seq, alphabet, code, used = encode_state_matrix(
        X,
        nr=options.get("nr", "*"),
        void=options.get("void", "%"),
    )
    sts_options = {key: value for key, value in options.items() if key in SEQ_KWARGS}
    for name in ("nr", "void"):
        text = used[name]
        if text is not None:
            sts_options[name] = text
        sts_options[f"{name}_code"] = None if text is None else code[text]
    sts_options["labels"], sts_options["cpal"] = seqdef_labels_and_palette(
        alphabet,
        nr=used["nr"],
        void=used["void"],
        labels=options.get("labels"),
        cpal=options.get("cpal"),
        missing_color=options.get("missing_color", "darkgrey"),
    )
    return STSList(data=seq, alphabet=alphabet, **sts_options)


def seqdef_labels_and_palette(
    alphabet: list[str],
    *,
    nr: str | None,
    void: str | None,
    labels: str | Sequence[object] | None = None,
    cpal: str | Sequence[object] | None = None,
    missing_color: str = "darkgrey",
) -> tuple[list[str], list[str]]:
    """One label and one colour per codebook entry, paired as ``seqdef`` pairs them.

    ``alphabet`` lists the real states, then the markers the data use (``nr``,
    ``void``; ``None`` when unused). As R's ``seqdef`` (and :func:`traminex.seqdef`),
    ``labels`` and ``cpal`` give one entry per real state, or one per codebook
    entry; the defaults are the states themselves and ``seqdef``'s palette, and a
    marker takes its own text as label and ``missing_color`` (nr) or white (void)
    as colour (TX-AUD-20261019-035, TX-REP-20261019-013).
    """
    markers = [text for text in (nr, void) if text is not None]
    real = [state for state in alphabet if state not in markers]
    marker_labels = {text: text for text in markers}
    marker_colours = {text: missing_color if text == nr else "white" for text in markers}

    def paired(
        values: str | Sequence[object] | None,
        defaults: list[str],
        extra: dict[str, str],
        message: str,
    ) -> list[str]:
        if values is None:
            own = dict(zip(real, defaults))
        else:
            listed = [str(v) for v in ([values] if isinstance(values, str) else values)]
            if len(listed) == len(alphabet):
                return listed
            if len(listed) != len(real):
                raise ValueError(message)
            own = dict(zip(real, listed))
        return [own[state] if state in own else extra[state] for state in alphabet]

    return (
        paired(
            labels,
            list(real),
            marker_labels,
            "\n [!] number of labels must equal number of states in the alphabet",
        ),
        paired(
            cpal,
            _default_seqdef_palette(len(real)),
            marker_colours,
            "\n [!] number of colors in 'cpal' must equal length of alphabet",
        ),
    )
