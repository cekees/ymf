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

# ---------------------------------------------------------------------------
# Shared enums (unchanged from v0.1)
# ---------------------------------------------------------------------------

ProvenanceEnum = Enum(
    ["llm_derived", "human_specified", "proteus_derived", "human_edited"]
)

FEFamilyEnum = Enum(["CG", "DG", "RT", "BDM", "RTF", "N1curl"])

BCTypeEnum = Enum(["dirichlet", "neumann", "robin", "periodic", "hydrograph"])

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
        Optional("tolerance"): Float(),
        Optional("max_iterations"): Int(),
    }
)

AnalyticalEntry = Map(
    {
        "name": Str(),
        "provenance": ProvenanceEnum,
        Optional("from_weak_form"): Any(),
        "method": Str(),
        "solution": Map(
            {
                "formula": Str(),
                Optional("notes"): Str(),
                Optional("derivation_source"): Str(),
                Optional("domain_restriction"): Str(),
                # Machine-readable: unknown name -> expression in x, y, z, t
                # and the coefficients (a vector unknown takes "(a, b)").
                # ymf.symbolic evaluates these; ``formula`` stays for people.
                Optional("expressions"): MapPattern(Str(), Str()),
            }
        ),
    }
)

DiscretizationEntry = Map(
    {
        "name": Str(),
        "provenance": ProvenanceEnum,
        Optional("from_weak_form"): Any(),
        Optional("finite_element"): FiniteElementDef,
        Optional("discretization"): Any(),
        "solver": SolverConfig,
        Optional("flux_method"): Str(),
        Optional("numerical_flux"): Str(),
        Optional("discretization_details"): Str(),
        Optional("proteus_modules"): Map({"physics": Str(), "numerics": Str()}),
        # --- v0.2 additions ---
        Optional("bmi_interface"): Bool(),
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
        Optional("units"): Str(),
        Optional("source"): Str(),
        Optional("source_url"): Str(),
        Optional("temporal"): Str(),
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

CharacteristicScalesDef = Map(
    {
        Optional("length"): ScaleDef,
        Optional("velocity"): ScaleDef,
        Optional("time"): ScaleDef,
        Optional("density"): ScaleDef,
        Optional("viscosity"): ScaleDef,
        Optional("temperature"): ScaleDef,
    }
)

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
        "strong_form_expression": Str(),
        "domain": Str(),
        "boundary_regions": Seq(Map({
            "name": Str(),
            "geometry": Str(),
            # Machine-readable membership test, e.g. "x = 0 or x = 4",
            # "x = 0 and y = 0"; ``geometry`` stays for people.
            Optional("where"): Str(),
        })),
        Optional("initial_conditions"): Seq(InitialConditionDef),
        Optional("boundary_conditions"): Seq(BoundaryConditionDef),
        Optional("known_analytical_solution"): Map(
            {"formula": Str(), "reference": Str()}
        ),
        Optional("coefficients"): MapPattern(Str(), Any()),
        Optional("dimensional_check"): DimensionalCheckDef,
        # Machine-readable equations, one string per (vector or scalar)
        # equation, in the operator syntax of ymf.symbolic:
        # "rho*dt(v) + div(rho*outer(v, v)) - div(mu*grad(v)) + grad(p) = f".
        # Ordered like the unknowns; strong_form_expression stays for people.
        Optional("equations"): Seq(Str()),
    }
)

# ---------------------------------------------------------------------------
# Problem (updated: + characteristic_scales, units, unit_systems, mesh_generation)
# ---------------------------------------------------------------------------

ProblemDef = Map(
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
        # Machine-readable domain. For now an axis-aligned box, whose
        # length also fixes the space dimension.
        Optional("geometry"): Map({
            "box": Map({"lower": Seq(Float()), "upper": Seq(Float())}),
        }),
    }
)

# ---------------------------------------------------------------------------
# Top-level YMF schema
# ---------------------------------------------------------------------------

YMF_SCHEMA = Map(
    {
        "Problem": ProblemDef,
        "solution_paths": DiscretizationCategory,
        Optional("vvuq"): Map(
            {Optional("verification"): Any(), Optional("validation"): Any()}
        ),
        Optional("archive"): Any(),
        Optional("Xdmf"): Any(),
        # File-level keys, see ymf.compose. "extends" names the document
        # this one builds on; "kind: model" marks a file that states the
        # physics but is not yet a well-posed problem; "composition" is
        # written by the loader into a composed document.
        Optional("extends"): Str(),
        Optional("kind"): Enum(["model", "problem"]),
        Optional("composition"): Any(),
    }
)


def _partial(validator):
    """The same validator with every mapping key made optional.

    Used for files that hold only part of a problem: a model, or a problem
    that extends one. Only mappings reached through mapping keys are
    relaxed. Entries of lists (an unknown, a weak form, a discretization)
    and of pattern maps keep their required keys, since an entry that is
    present at all should be complete.
    """
    if isinstance(validator, Map):
        return Map({
            (key if isinstance(key, Optional) else Optional(key)): _partial(value)
            for key, value in validator._validator.items()
        })
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
    a ``kind: model`` file is validated as a partial document; see
    :mod:`ymf.compose`. Raises ``strictyaml.YAMLValidationError`` on schema
    violations, naming the file that holds the error.
    """
    from ymf.compose import load_composed  # ymf.compose imports this module

    return load_composed(path)


def validate_ymf(text: str) -> AnyType:
    """Validate a YMF document already loaded as a string; returns plain data.

    Does not follow ``extends``, which needs a file to resolve paths
    against; use :func:`load_ymf` for that.
    """
    parsed = strictyaml.dirty_load(
        text, YMF_SCHEMA, allow_flow_style=_ALLOW_FLOW_STYLE
    )
    return parsed.data
