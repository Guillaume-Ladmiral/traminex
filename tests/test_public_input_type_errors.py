"""A wrong-type input raises TypeError, not AttributeError (TX-AUD-20261010-029, TX-REP-20261018-001).

Called with a plain array, sixteen public callables reached for an attribute of
the object they expect (``_canonicalized``, ``channels``, ``id``, ``label``,
``columns``, ``model``, ``info``, ``adjustment``) and raised AttributeError,
where the rest raise TypeError naming the expected type: seqidegrad, seqibad and
seqipos (TX-AUD-20261010-029) and thirteen event, tree and multichannel
functions (TX-REP-20261018-001). The census calls every public callable whose
first parameter is positional with a plain 2-D integer array.
"""

from __future__ import annotations

import inspect
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

import traminex as tx  # noqa: E402
from traminex.api_manifest import PUBLIC_API  # noqa: E402


def _positional_callables() -> list[str]:
    names = []
    for name in PUBLIC_API:
        obj = getattr(tx, name)
        if inspect.isclass(obj) or not callable(obj):
            continue
        try:
            params = list(inspect.signature(obj).parameters.values())
        except (TypeError, ValueError):
            continue
        if params and params[0].kind in (
            params[0].POSITIONAL_ONLY,
            params[0].POSITIONAL_OR_KEYWORD,
        ):
            names.append(name)
    return names


@pytest.mark.parametrize("name", _positional_callables())
def test_a_plain_array_never_raises_attribute_error(name: str) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            getattr(tx, name)(np.array([[0, 0, 1], [1, 1, 0]]))
        except AttributeError as exc:  # the defect
            pytest.fail(f"{name} raised AttributeError: {exc}")
        except Exception:  # noqa: BLE001 - any other outcome is the function's own contract
            pass
        finally:
            plt.close("all")


#: The census above passes the array alone, so a callable with another required parameter
#: stopped at the missing argument before it read the array; with those arguments supplied
#: these sixteen still raised AttributeError (TX-AUD-20261019-036).
SUPPLIED = (
    "seqMD_gridsearch_weights",
    "seqMD_optimize_weights",
    "seqrecode",
    "seqimplic",
    "seqsha",
    "seqstart",
    "cluster_representatives",
    "cluster_summary",
    "frequent_event_pairs_timegap",
    "seqtree_cv_fit",
    "seqtree_fit",
    "seqtree_predict",
    "seqtree2dot",
    "disstree2dot_explore",
    "disstree2dotp",
    "disstree_predict",
)


def _other_required_arguments(name: str) -> tuple[list, dict]:
    """A stand-in for each further required parameter: 1, or [0, 1] for a group (the audit's probe)."""
    positional: list = []
    keywords: dict = {}
    for p in list(inspect.signature(getattr(tx, name)).parameters.values())[1:]:
        if p.default is not inspect.Parameter.empty:
            continue
        value = [0, 1] if p.name in ("groups", "group") else 1
        if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD):
            positional.append(value)
        elif p.kind == p.KEYWORD_ONLY:
            keywords[p.name] = value
    return positional, keywords


def _call_with_a_plain_array(name: str) -> None:
    positional, keywords = _other_required_arguments(name)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            getattr(tx, name)(np.array([[0, 0, 1], [1, 1, 0]]), *positional, **keywords)
        finally:
            plt.close("all")


@pytest.mark.parametrize(
    "name", [name for name in _positional_callables() if any(_other_required_arguments(name))]
)
def test_a_plain_array_with_the_other_arguments_never_raises_attribute_error(name: str) -> None:
    try:
        _call_with_a_plain_array(name)
    except AttributeError as exc:  # the defect
        pytest.fail(f"{name} raised AttributeError: {exc}")
    except Exception:  # noqa: BLE001 - any other outcome is the function's own contract
        pass


@pytest.mark.parametrize("name", SUPPLIED)
def test_the_object_argument_is_checked_first(name: str) -> None:
    with pytest.raises(TypeError, match="must be"):
        _call_with_a_plain_array(name)
