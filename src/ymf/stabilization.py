"""What a ``stabilization`` entry says: the method, its length scale, and the shock-capturing viscosity with its parameters.

A sub-problem (or a split discretization's entry, as the default of its sub-problems) states its stabilization as a word
or as a block::

    stabilization: SUPG                                  # the method alone
    stabilization:
      method: SUPG
      tau: metric                    # the length scale: the element metric G (proteus useMetrics 1) or the size h (0)
      lagged: false                  # the subgrid (tau) error taken from the previous step (SubgridError lag)
      shock_capturing:
        factor: 0.5                  # the shock-capturing viscosity's coefficient (shockCapturingFactor, sc_c)
        beta: 1.5                    # the exponent of the gradient it is scaled by (sc_beta; 1 = classical)
        reference: 1.0               # the reference value of the unknown (sc_uref)
        lagged: true                 # the viscosity taken from the previous step (ShockCapturing lag)
    stabilization: {method: SUPG, shock_capturing: none}     # none, stated: the solver adds no shock capturing
    stabilization: none                                  # nothing: no method, no shock capturing (MCorr, a projection)
    stabilization: {shock_capturing: {factor: 0.9}}      # shock capturing with no method (proteus's RDLS block of 0.2)

What it is for. The generator (C19) turned VOF and NCLS into a functor each, but the SUPG tau, the shock capturing and
their lag are the engine template's, and the template took them from the model's name (``dambreak.ymf`` stated none: the
generator spike's ``STAB`` was declared in ``make_inputs.py``). A spec that states them lets the engine read them.

The words are checked on load (:func:`check`, called by ``ymf.vocabulary``); :func:`resolve` normalizes a block (the
alias ``shock_capturing`` of the method word becomes ``residual_shock_capturing``, and the keys proteus's dispatcher
already reads become ymf's: ``lag`` -> ``lagged``, ``sc_beta`` -> ``beta``, ``sc_uref`` -> ``reference``; the bare word
becomes a block) and adds nothing: **a key the spec does not state stays unstated**, because the value a solver uses then is its own default
(``proteus.DispatchTable``'s cells), and two defaults, one in ymf and one in the solver, would be two meanings.
:func:`missing` lists what a reader that must not guess (the generator) still needs stated.

Which keys a method takes:

* ``none``: no residual-based method. It takes only ``shock_capturing``: ``stabilization: none`` (or ``method: none``)
  says that the sub-problem is not stabilized at all, shock capturing included; ``{method: none, shock_capturing:
  {factor: 0.9}}`` says shock capturing and nothing else. A block that states ``shock_capturing`` and no ``method`` is
  exactly that second one (ymf 0.2 read such a block, proteus's RDLS stated one; 0.3.0 refused it; resolve gives it
  ``method: none``). This is a meaning, not a default: ymf never picks a method for a block that does not name one.
* ``SUPG``, ``ASGS``: ``tau``, ``lagged``, ``shock_capturing`` (the residual-based viscosity added to the method).
* ``residual_shock_capturing``: the shock-capturing viscosity alone: ``tau``, ``shock_capturing`` (its block says how much;
  unstated, the solver's own coefficient applies).
* ``entropy_viscosity``, ``graph_viscosity_smoothness``, ``graph_viscosity_entropy_residual``: no keys here. Their
  parameters (the entropy-viscosity coefficient, the lumped mass) are the scheme's options in the ``discretization:`` block.

Stdlib only, like the rest of the spec core.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from ymf.vocabulary import SHOCK_CAPTURING_KEYS, TAU_FORMS, VocabularyError, canonical

__all__ = ["StabilizationError", "TAU_FORMS", "TAKES", "SHOCK_CAPTURING_KEYS", "resolve", "check", "for_subproblem", "missing",
           "generator_needs"]


class StabilizationError(ValueError):
    """A stabilization a spec states that does not make sense."""


#: method -> the keys of the block it reads (beside ``method``)
TAKES: Dict[str, tuple] = {
    "none": ("shock_capturing",),
    "SUPG": ("tau", "lagged", "shock_capturing"),
    "ASGS": ("tau", "lagged", "shock_capturing"),
    "residual_shock_capturing": ("tau", "shock_capturing"),
    "entropy_viscosity": (),
    "graph_viscosity_smoothness": (),
    "graph_viscosity_entropy_residual": (),
}

def _number(value: Any, where: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        raise StabilizationError("%s: %r is not a number" % (where, value)) from None


def _bool(value: Any, where: str) -> bool:
    if isinstance(value, bool):
        return value
    if str(value).lower() in ("true", "1"):
        return True
    if str(value).lower() in ("false", "0"):
        return False
    raise StabilizationError("%s: %r is not true or false" % (where, value))


def _word(kind: str, word: Any, where: str) -> str:
    try:
        return canonical(kind, word)
    except VocabularyError as exc:
        raise StabilizationError("%s: %s" % (where, exc)) from None


def _method(word: Any, where: str) -> str:
    return _word("stabilization", word, where)


def _keys(block: Dict[str, Any], kind: str, where: str) -> Dict[str, Any]:
    """The block with its keys in ymf's words; a word and its alias together are refused."""
    out: Dict[str, Any] = {}
    for key, value in block.items():
        try:
            word = canonical(kind, key)
        except VocabularyError:
            word = key                      # not a word of the kind: the caller names what it reads
        if word in out:
            raise StabilizationError("%s: %s is stated twice (%s is the old spelling of %s)"
                                     % (where, word, key if key != word else "an alias", word))
        out[word] = value
    return out


def _shock_capturing(block: Any, where: str) -> Optional[Dict[str, Any]]:
    if isinstance(block, str):
        if block != "none":
            raise StabilizationError("%s: %r is not a block; say none, or give factor, beta, reference, lagged"
                                     % (where, block))
        return None
    if not isinstance(block, dict):
        raise StabilizationError("%s: a block of factor, beta, reference, lagged, or none" % (where,))
    block = _keys(block, "shock_capturing_key", where)
    unknown = sorted(set(block) - set(SHOCK_CAPTURING_KEYS))
    if unknown:
        raise StabilizationError("%s: no key %s in shock_capturing; its keys are %s"
                                 % (where, ", ".join(unknown), ", ".join(SHOCK_CAPTURING_KEYS)))
    if "factor" not in block:
        raise StabilizationError("%s: factor is required (a block that states none says none)" % (where,))
    out: Dict[str, Any] = {}
    out["factor"] = _number(block["factor"], where + ".factor")
    if out["factor"] < 0:
        raise StabilizationError("%s.factor: %g is negative" % (where, out["factor"]))
    if "beta" in block:
        out["beta"] = _number(block["beta"], where + ".beta")
        if not 1.0 <= out["beta"] <= 2.0:
            raise StabilizationError("%s.beta: %g is outside [1, 2] (1 is the classical scaling by |grad u|, "
                                     "2 makes the viscosity independent of it)" % (where, out["beta"]))
    if "reference" in block:
        out["reference"] = _number(block["reference"], where + ".reference")
        if out["reference"] <= 0:
            raise StabilizationError("%s.reference: %g is not positive" % (where, out["reference"]))
    if "lagged" in block:
        out["lagged"] = _bool(block["lagged"], where + ".lagged")
    return out


def resolve(block: Any, where: str = "stabilization") -> Dict[str, Any]:
    """The stabilization as a block with ymf's words and every stated key kept: ``{method, tau?, lagged?,
    shock_capturing?}``, ``shock_capturing`` a dict, or ``None`` when the spec says ``none``. Nothing is added, with one
    reading of what is stated: a block that gives ``shock_capturing`` and no ``method`` has ``method: none`` (shock
    capturing and no residual-based method), and the bare word ``none`` is ``{"method": "none"}``, which under
    :data:`TAKES` ``none`` states no shock capturing either. Raises :class:`StabilizationError`."""
    if isinstance(block, str):
        block = {"method": block}
    if not isinstance(block, dict):
        raise StabilizationError("%s: say the method (stabilization: SUPG, or none), or give a block with method:"
                                 % (where,))
    block = _keys(block, "stabilization_key", where)
    if "method" not in block:
        if "shock_capturing" not in block:
            raise StabilizationError(
                "%s: say the method (stabilization: SUPG, or none), or give a block with method: (a block with only "
                "shock_capturing has no method: none)" % (where,))
        others = sorted(set(block) - {"shock_capturing"})
        if others:
            raise StabilizationError(
                "%s: no method is stated, so the block is shock capturing and nothing else (method: none), and "
                "%s %s read by a method: say which (method: SUPG, ASGS or residual_shock_capturing)"
                % (where, ", ".join(others), "is" if len(others) == 1 else "are"))
        block = dict(block, method="none")
    method = _method(block["method"], where + ".method")
    takes = TAKES[method]
    unread = sorted(set(block) - {"method"} - set(takes))
    if unread:
        raise StabilizationError(
            "%s: %s does not read %s; %s" % (
                where, method, ", ".join(unread),
                ("its keys are method, " + ", ".join(takes)) if takes else
                "it has no keys here (its parameters are the scheme's options in the discretization: block)"))
    out: Dict[str, Any] = {"method": method}
    if "tau" in block:
        tau = str(block["tau"])
        if tau not in TAU_FORMS:
            raise StabilizationError("%s.tau: %s is not one of %s" % (where, tau, ", ".join(TAU_FORMS)))
        out["tau"] = tau
    if "lagged" in block:
        out["lagged"] = _bool(block["lagged"], where + ".lagged")
    if "shock_capturing" in block:
        out["shock_capturing"] = _shock_capturing(block["shock_capturing"], where + ".shock_capturing")
    return out


def check(block: Any, where: str = "stabilization") -> Dict[str, Any]:
    """:func:`resolve`, under the name vocabulary checks use."""
    return resolve(block, where)


def for_subproblem(entry: Dict[str, Any], sub: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """A sub-problem's resolved stabilization: its own block, or else its entry's (the default of every sub-problem
    of the entry); a sub-problem that states one replaces the entry's whole, key by key merging would make two
    blocks say a third thing. None when neither states one."""
    block = sub.get("stabilization")
    if block is None:
        block = entry.get("stabilization")
    if block is None:
        return None
    return resolve(block, "sub-problem %s stabilization" % (sub.get("name", "?"),))


def missing(block: Any) -> List[str]:
    """The keys a reader that must not fall back on a default still needs, as dotted names: what
    :func:`generator_needs` asks for, less what the block states. ``[]`` if it is complete."""
    r = resolve(block)
    return [need for need in generator_needs(r["method"], r) if not _has(r, need)]


def generator_needs(method: str, resolved: Optional[Dict[str, Any]] = None) -> List[str]:
    """What the kernel generator reads of a stabilization of ``method``: ``tau`` (and ``lagged`` for SUPG and ASGS),
    and the shock-capturing block's four keys (or its explicit ``none``). Method ``none`` needs no statement that there
    is no shock capturing (it says so); only a shock-capturing block it does state must be complete."""
    takes = TAKES[method]
    needs = [k for k in ("tau", "lagged") if k in takes]
    if "shock_capturing" in takes:
        sc = (resolved or {}).get("shock_capturing", "unstated")
        if method != "none":
            needs.append("shock_capturing")
        if isinstance(sc, dict):
            needs += ["shock_capturing." + k for k in SHOCK_CAPTURING_KEYS]
    return needs


def _has(resolved: Dict[str, Any], dotted: str) -> bool:
    here: Any = resolved
    for part in dotted.split("."):
        if not isinstance(here, dict) or part not in here:
            return False
        here = here[part]
    return True
