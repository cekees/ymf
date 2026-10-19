"""strictyaml schema validators for YMF v0.2.

This module implements the schema described in
``docs/ymf-schema-v0.2-delta.md``, merged on top of the v0.1 schema in
``docs/ymf-schema.md``. Relative to v0.1, ``unknowns`` widens its
accepted type (bare string OR structured map), and the strong form's
``unknown_provenance`` is replaced by ``provenance``, the key every other
block uses; a structured unknown may carry its own ``provenance`` to
override it. Every other v0.2 addition is optional.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any as AnyType

import strictyaml
from strictyaml import (
    Any,
    Bool,
    Enum,
    Float,
    Int,
    Map,
    MapPattern,
    Optional,
    Seq,
    Str,
)

from ymf.discretizations import METHODS
from ymf.state import STATES

# ---------------------------------------------------------------------------
# Shared enums (unchanged from v0.1)
# ---------------------------------------------------------------------------

ProvenanceEnum = Enum(
    ["llm_derived", "human_specified", "proteus_derived", "human_edited"]
)

FEFamilyEnum = Enum(["CG", "DG", "RT", "BDM", "RTF", "N1curl"])

# The mathematical kind only. How the value is *supplied* -- a number, a
# formula, or a table -- is orthogonal, and lives in value/formula/series
# below. "hydrograph" used to be listed here, which tangled the two: a
# discharge hydrograph is a tabulated neumann condition, and a stage
# hydrograph a tabulated dirichlet one, which naming it as a type made
# impossible to write.
BCTypeEnum = Enum(["dirichlet", "neumann", "robin", "periodic"])

ICTypeEnum = Enum(["constant", "function", "data_source"])

UnitConventionEnum = Enum(["SI", "CGS", "US", "dimensionless"])

DimCheckMethodEnum = Enum(["automatic", "manual", "off"])

# ---------------------------------------------------------------------------
# v0.1 building blocks (unchanged)
# ---------------------------------------------------------------------------

WeakFormDef = Map(
    {
        "label": Str(),
        "provenance": ProvenanceEnum,
        "derivation": Str(),
        "solution_spaces": Map({"trial": Str(), "test": Str()}),
        Optional("inf_sup_stability"): Str(),
        Optional("existence_note"): Str(),
        Optional("uniqueness_note"): Str(),
        Optional("regularity_note"): Str(),
        "bilinear": Str(),
        "linear": Str(),
        Optional("boundary_conditions"): Str(),
        Optional("stabilization_method"): Str(),
        Optional("notes"): Str(),
    }
)

FiniteElementDef = Map(
    {
        Optional("velocity"): Map({"family": FEFamilyEnum, "order": Int()}),
        Optional("pressure"): Map({"family": FEFamilyEnum, "order": Int()}),
        Optional("fields"): Map({"family": FEFamilyEnum, "order": Int()}),
    }
)

SolverConfig = Map(
    {
        Optional("type"): Enum(
            ["linear", "nonlinear", "time_marching", "saddle_point"]
        ),
        Optional("nonlinear_solver"): Enum(["newton", "broyden", "line_search"]),
        Optional("linear_solver"): Enum(["petsc", "mumps", "umfpack", "superlu"]),
        # PETSc options, without the leading dash, for linear_solver petsc
        # or mumps: {ksp_type: cg, pc_type: gamg}. Values are passed as
        # written; an option that takes none is given "".
        Optional("petsc_options"): MapPattern(Str(), Str()),
        Optional("tolerance"): Float(),
        Optional("max_iterations"): Int(),
    }
)

# No from_weak_form here. An analytical solution solves the strong form;
# the weak formulation is the business of theory and numerics. (Non-classical
# solutions built directly from the integral form of a conservation law are a
# real thing, but they are not what this key meant, and nothing used it.)
AnalyticalEntry = Map(
    {
        "name": Str(),
        "provenance": ProvenanceEnum,
        "method": Str(),
        "solution": Map(
            {
                "formula": Str(),
                Optional("notes"): Str(),
                Optional("derivation_source"): Str(),
                Optional("domain_restriction"): Str(),
            }
        ),
    }
)

DiscretizationEntry = Map(
    {
        "name": Str(),
        "provenance": ProvenanceEnum,
        # A weak form's label. Integer indices into weak_forms used to be
        # accepted too; they are not, because "remove:" can delete a weak
        # form and silently renumber every index after it. May be omitted
        # when the document has exactly one weak form -- ymf.references
        # fills it in during composition.
        Optional("from_weak_form"): Str(),
        Optional("finite_element"): FiniteElementDef,
        Optional("discretization"): Any(),
        "solver": SolverConfig,
        Optional("flux_method"): Str(),
        Optional("numerical_flux"): Str(),
        Optional("discretization_details"): Str(),
        Optional("proteus_modules"): Map({"physics": Str(), "numerics": Str()}),
        # --- v0.2 additions ---
        Optional("bmi_interface"): Bool(),
        # How the model is reduced to a computable problem, and so which of
        # the blocks below this entry needs: mesh (the default), neural, or
        # particle. See ymf.discretizations, which also says what each
        # method still lacks; the fields here are optional so that a
        # missing one is reported with a hint, not as a schema error.
        Optional("method"): Enum(list(METHODS)),
        # The mesh a run uses: cells along each side of the domain's box,
        # or h, the element size of an unstructured mesh of any other
        # region (one of the two), and levels of uniform refinement of that
        # coarsest mesh. A list is a study: each value is one realization
        # (see ymf.closure).
        Optional("mesh"): Map({
            Optional("cells"): Int() | Seq(Int()),
            Optional("h"): Float() | Seq(Float()),
            Optional("levels"): Int() | Seq(Int()),
        }),
        # method: neural -- a network trained on the equations (PINN, weak
        # adversarial network). A schema so far: no outputs are planned.
        Optional("neural"): Map({
            Optional("width"): Int() | Seq(Int()),
            Optional("depth"): Int() | Seq(Int()),
            Optional("points"): Int() | Seq(Int()),
            Optional("activation"): Str(),
            Optional("training"): Any(),
        }),
        # method: particle -- a meshless particle method (RKPM, SPH). A
        # schema so far: no outputs are planned.
        Optional("particle"): Map({
            Optional("spacing"): Float() | Seq(Float()),
            Optional("kernel"): Str(),
            Optional("support"): Float(),
        }),
        Optional("roughness_source"): Map(
            {
                "type": Str(),
                Optional("land_cover_data"): Map(
                    {
                        "source": Str(),
                        Optional("resolution"): Str(),
                        Optional("url"): Str(),
                    }
                ),
                Optional("reference"): Str(),
                Optional("lookup_table"): Str(),
            }
        ),
    }
)

DiscretizationCategory = Map(
    {
        "analytical": Seq(AnalyticalEntry),
        "discretizations": Seq(DiscretizationEntry),
    }
)

PhysicalModelDef = Map(
    {
        "provenance": ProvenanceEnum,
        Optional("processes"): Seq(Str()),
        Optional("assumptions"): Seq(Str()),
        Optional("source_document"): Str(),
        Optional("notes"): Str(),
    }
)

# ---------------------------------------------------------------------------
# v0.2 additions
# ---------------------------------------------------------------------------

# --- unknowns: BREAKING widen from Seq(Str()) to Seq(Str() | Map(...)) ---
UnknownDef = Str() | Map(
    {
        "name": Str(),
        Optional("units"): Str(),
        Optional("std_name"): Str(),
        # Overrides strong_form.provenance for this unknown alone.
        Optional("provenance"): ProvenanceEnum,
        # 0 for a scalar (the default), 1 for a vector with one component
        # per space dimension.
        Optional("rank"): Int(),
    }
)

InitialConditionDef = Map(
    {
        "field": Str(),
        "type": ICTypeEnum,
        Optional("value"): Float(),
        Optional("formula"): Str(),
        Optional("source"): Str(),
        Optional("source_url"): Str(),
        Optional("temporal"): Str(),
        Optional("units"): Str(),
    }
)

BoundaryConditionDef = Map(
    {
        "region": Str(),
        "variable": Str(),
        "type": BCTypeEnum,
        Optional("value"): Float(),
        # A value that varies along the boundary or in time: an expression
        # in x, y, z, t and the coefficients; "(a, b)" for a vector variable.
        Optional("formula"): Str(),
        # The value as tabulated data rather than a formula -- a measured or
        # design record, which is the one thing that cannot be written as an
        # expression. Pairs of [t, value], t ascending. Linear between
        # samples; outside the range the end value is held, because
        # extrapolating a record past its end is not something to do
        # quietly.
        #
        # Available to any type: with neumann this is a discharge
        # hydrograph, with dirichlet a stage hydrograph or a measured inlet
        # temperature.
        Optional("series"): Seq(Seq(Float())),
        # A robin condition, alpha*u + beta*(q.n) = g: alpha and beta here,
        # g in value or formula like any other right-hand side. Expressions,
        # so they may vary along the boundary.
        Optional("alpha"): Str(),
        Optional("beta"): Str(),
        Optional("units"): Str(),
        Optional("source"): Str(),
        Optional("source_url"): Str(),
        Optional("temporal"): Str(),
    }
)

# A domain is the notation short form -- "Ω = [0, 1] × [0, 1],  t ∈ (0, 50]",
# which means a box and is parsed by ymf.symbolic -- or a mapping tagged by
# kind. See ymf.domain, which checks the mapping against its kind; the schema
# only says what keys may appear, because strictyaml cannot express "these
# keys given that kind" and a message naming the kind is worth more than one
# naming the schema.
DomainDef = Str() | Map(
    {
        "kind": Enum(["box", "piecewise_linear", "discrete", "file"]),
        Optional("name"): Str(),
        Optional("dimension"): Int(),
        Optional("units"): Str(),
        # the time interval of a transient model, as the notation
        # form states it ("..., t ∈ (0, 50]")
        Optional("time"): Seq(Float()),
        # A name for each boundary, so a condition can say "bottom" instead
        # of 3. The integer is the flag Proteus and the mesh generators use.
        Optional("boundary_tags"): MapPattern(Str(), Int()),
        # kind: box
        Optional("lower"): Seq(Float()),
        Optional("upper"): Seq(Float()),
        # kind: piecewise_linear -- a boundary representation
        Optional("vertices"): Seq(Seq(Float())),
        Optional("segments"): Seq(Seq(Int())),
        Optional("facets"): Any(),
        Optional("volumes"): Any(),
        Optional("holes"): Seq(Seq(Float())),
        # A region carries its own flag; vertex/segment/facet flags stay
        # parallel arrays, which is what they are and what Proteus holds.
        Optional("regions"): Seq(Map({"point": Seq(Float()),
                                      Optional("flag"): Int()})),
        Optional("vertex_flags"): Seq(Int()),
        Optional("segment_flags"): Seq(Int()),
        Optional("facet_flags"): Seq(Int()),
        # kind: discrete -- already a mesh, in the archive's vocabulary
        Optional("topology"): Any(),
        Optional("geometry"): Any(),
        # kind: file -- someone else's language
        Optional("format"): Str(),
        Optional("path"): Str(),
    }
)

ScaleDef = Map(
    {
        Optional("value"): Float(),
        Optional("units"): Str(),
        Optional("derived"): Bool(),
        Optional("formula"): Str(),
    }
)

# A coefficient is a number, a map, or an expression -- tried in that order.
#
# Typing this mattered more than it looks. It used to be Any(), under which
# strictyaml returns every scalar as a string: "kappa: 1.0e-3" loaded as the
# string '1.0e-3'. Only compute_dimensionless_numbers coerced it, so every
# other consumer saw a string, and a formula coefficient that referred to a
# numeric one silently dropped out of scope, taking any number built on it.
#
# Float() first, so a numeric scalar is a number however it was written --
# and so strictyaml agrees with the pyyaml that reads a spec carried inside
# an archive. Str() last, as the fallback for an expression ("π²",
# "ρ/(2μ) - sqrt(ρ²/(4μ²) + 4π²)"), which is a coefficient too.
#
# The map form keeps every key optional on purpose: a model file declares
# "ρ: {units: kg/m3}" with no value, leaving the value to the problem that
# extends it.
CoefficientMapDef = Map(
    {
        Optional("value"): Float(),
        Optional("units"): Str(),
        Optional("formula"): Str(),
    }
)

CoefficientDef = Float() | CoefficientMapDef | Str()

# Any named scale, not a fixed set.
#
# This was a closed Map of six: length, velocity, time, density, viscosity,
# temperature. A problem with a diffusivity or a pressure scale had nowhere
# to put it and had to leave it out of its dimensionless numbers.
#
# Those six keep their meaning: characteristic_scale_values aliases them to
# the bare names L, U, rho, mu, T that the report's formulas use. Any other
# name resolves under itself, so "diffusivity: {value: 1.0e-3}" is usable in
# a formula as "diffusivity".
CharacteristicScalesDef = MapPattern(Str(), ScaleDef)

DimensionlessNumberDef = Map({"formula": Str(), Optional("value"): Float()})

DimensionalCheckDef = Map(
    {
        Optional("enabled"): Bool(),
        Optional("method"): DimCheckMethodEnum,
        Optional("tolerance"): Float(),
    }
)

UnitSystemDef = Map(
    {
        Optional("length"): Str(),
        Optional("mass"): Str(),
        Optional("time"): Str(),
        Optional("temperature"): Str(),
    }
)

UnitsDef = Map(
    {
        Optional("convention"): UnitConventionEnum,
        Optional("default"): Bool(),
    }
)

MeshGenerationDef = Map(
    {
        "source": Str(),
        Optional("source_url"): Str(),
        Optional("format"): Str(),
        Optional("vertical_datums"): Map(
            {Optional("horizontal"): Str(), Optional("vertical"): Str()}
        ),
        Optional("resolution_target"): Float(),
        Optional("filter"): Map(
            {
                Optional("method"): Str(),
                Optional("kernel"): Str(),
                Optional("sigma"): Float(),
            }
        ),
    }
)

# ---------------------------------------------------------------------------
# StrongFormDef (updated: widened unknowns + new IC/BC/dimensional_check)
# ---------------------------------------------------------------------------

StrongFormDef = Map(
    {
        # Covers the equations and the choice of unknowns. Earlier drafts
        # called this "unknown_provenance", which read as "provenance unknown".
        "provenance": ProvenanceEnum,
        "unknowns": Seq(UnknownDef),  # CHANGED in v0.2: was Seq(Str())
        "equation_formulation": Str(),
        # Superseded by ``equations``, which is written in the same notation
        # and read by ymf.symbolic; kept so older documents still validate.
        Optional("strong_form_expression"): Str(),
        # "Ω = [0, 1] × [0, 1]", optionally ", t ∈ (0, 50]" for the time
        # interval of a transient problem.
        "domain": DomainDef,
        # geometry is a membership test in the spec notation: "∂Ω" (all of
        # it), "y = 0 or y = H", "x = 0 and y = 0".
        # "where" is a membership predicate, not geometry: "∂Ω" (all of
        # it), "y = 0 or y = H", "x = 0 and y = 0". adr_problem already
        # calls this thing "where" in the conditions it returns.
        "boundary_regions": Seq(Map({"name": Str(), "where": Str()})),
        Optional("initial_conditions"): Seq(InitialConditionDef),
        Optional("boundary_conditions"): Seq(BoundaryConditionDef),
        Optional("known_analytical_solution"): Map(
            {"formula": Str(), "reference": Str()}
        ),
        Optional("coefficients"): MapPattern(Str(), CoefficientDef),
        Optional("dimensional_check"): DimensionalCheckDef,
        # The equations, one string per (vector or scalar) equation, in the
        # notation of ymf.symbolic.notation -- "∂T/∂t = ∇·(κ∇T)  in Ω",
        # "ρ ∂v/∂t + ∇·(ρ v⊗v) − μΔv + ∇p = f" -- ordered like the unknowns.
        Optional("equations"): Seq(Str()),
    }
)

# ---------------------------------------------------------------------------
# Problem (updated: + characteristic_scales, units, unit_systems, mesh_generation)
# ---------------------------------------------------------------------------

ModelDef = Map(
    {
        "name": Str(),
        Optional("physical_model"): PhysicalModelDef,
        Optional("strong_form"): StrongFormDef,
        Optional("weak_forms"): Seq(WeakFormDef),
        Optional("characteristic_scales"): CharacteristicScalesDef,
        Optional("dimensionless_numbers"): MapPattern(Str(), DimensionlessNumberDef),
        Optional("units"): UnitsDef,
        Optional("unit_systems"): MapPattern(Str(), UnitSystemDef),
        Optional("mesh_generation"): MeshGenerationDef,
    }
)

# ---------------------------------------------------------------------------
# Top-level YMF schema
# ---------------------------------------------------------------------------

YMF_SCHEMA = Map(
    {
        "model": ModelDef,
        "solution_paths": DiscretizationCategory,
        Optional("vvuq"): Map(
            {Optional("verification"): Any(), Optional("validation"): Any()}
        ),
        Optional("archive"): Any(),
        Optional("xdmf"): Any(),
        # File-level keys, see ymf.compose. "extends" names the document
        # this one builds on; "ymf_state" says how far the document has got
        # (abstract, well-posed, discretized, closed: see ymf.state), and is
        # checked against its content; "composition" is written by the
        # loader into a composed document.
        Optional("extends"): Str(),
        Optional("ymf_state"): Enum(list(STATES)),
        # "generated" says a tool wrote this file, from what, and digests
        # of what it read and what it wrote: see ymf.generated.
        Optional("generated"): Map({
            "by": Str(),
            "reads": Enum(["model_without_weak_forms", "model_and_analytical"]),
            "from": Str(),
            "content": Str(),
        }),
        Optional("composition"): Any(),
        # "remove" lists paths this file drops from what it extends --
        # "Problem.physical_model.assumptions[incompressible]". Composition
        # could only add and change before, so dropping an assumption or a
        # field took a new model file even when the change was one line.
        Optional("remove"): Seq(Str()),
    }
)


def _partial(validator, keep=()):
    """The same validator with every mapping key made optional.

    Used for files that hold only part of a problem: a model, or a problem
    that extends one. Entries of lists are relaxed too, except for the key
    that names them (``name`` or ``label``): a file may change one field of
    an inherited discretization by naming it and giving only that field.
    The composed document is validated in full, so an entry that is new
    and incomplete is still caught there.
    """
    if isinstance(validator, Map):
        return Map({
            (key if isinstance(key, Optional) or key in keep else Optional(key)):
                _partial(value)
            for key, value in validator._validator.items()
        })
    if isinstance(validator, Seq) and isinstance(validator._validator, Map):
        return Seq(_partial(validator._validator, keep=("name", "label")))
    return validator


#: :data:`YMF_SCHEMA` with every section and field optional, for model files
#: and for the individual files of an ``extends`` chain.
PARTIAL_YMF_SCHEMA = _partial(YMF_SCHEMA)


# YMF documents and examples use YAML flow-style collections extensively
# (e.g. ``unknowns: [u, p]``, ``fields: {family: CG, order: 2}``).
# strictyaml disallows flow style by default (``allow_flow_style=False``)
# as a style-consistency feature; YMF explicitly opts back in since flow
# style is idiomatic for short inline lists/maps throughout the schema.
_ALLOW_FLOW_STYLE = True


def load_ymf(path: str | Path) -> AnyType:
    """Load and validate a YMF YAML file, returning a plain Python object.

    A file with ``extends`` is composed with the documents it builds on, and
    a ``ymf_state: abstract`` file is validated as a partial document; see
    :mod:`ymf.compose` and :mod:`ymf.state`. Raises ``strictyaml.YAMLValidationError`` on schema
    violations, naming the file that holds the error.
    """
    from ymf.compose import load_composed  # ymf.compose imports this module

    return load_composed(path)


def validate_ymf(text: str) -> AnyType:
    """Validate a YMF document already loaded as a string; returns plain data.

    Does not follow ``extends``, which needs a file to resolve paths
    against; use :func:`load_ymf` for that.
    """
    from ymf.compose import _KEY_ALIASES

    import re as _re
    for was, now in _KEY_ALIASES.items():
        text = _re.sub(r"^%s:" % was, "%s:" % now, text, flags=_re.M)
    parsed = strictyaml.dirty_load(
        text, YMF_SCHEMA, allow_flow_style=_ALLOW_FLOW_STYLE
    )
    return parsed.data
