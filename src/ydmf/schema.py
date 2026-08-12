"""strictyaml schema validators for YDMF v0.2.

This module implements the schema described in
``docs/ydmf-schema-v0.2-delta.md``, merged on top of the v0.1 schema in
``docs/ydmf-schema.md``. Only ``unknowns`` widens its accepted type
(bare string OR structured map) relative to v0.1; every other v0.2
addition is optional and preserves v0.1 documents unchanged.
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
        "unknowns": Seq(UnknownDef),  # CHANGED in v0.2: was Seq(Str())
        "unknown_provenance": ProvenanceEnum,
        "equation_formulation": Str(),
        "strong_form_expression": Str(),
        "domain": Str(),
        "boundary_regions": Seq(Map({"name": Str(), "geometry": Str()})),
        Optional("initial_conditions"): Seq(InitialConditionDef),
        Optional("boundary_conditions"): Seq(BoundaryConditionDef),
        Optional("known_analytical_solution"): Map(
            {"formula": Str(), "reference": Str()}
        ),
        Optional("coefficients"): MapPattern(Str(), Any()),
        Optional("dimensional_check"): DimensionalCheckDef,
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
    }
)

# ---------------------------------------------------------------------------
# Top-level YDMF schema
# ---------------------------------------------------------------------------

YDMF_SCHEMA = Map(
    {
        "Problem": ProblemDef,
        "solution_paths": DiscretizationCategory,
        Optional("vvuq"): Map(
            {Optional("verification"): Any(), Optional("validation"): Any()}
        ),
        Optional("archive"): Any(),
        Optional("Xdmf"): Any(),
    }
)


def load_ydmf(path: str | Path) -> AnyType:
    """Load and validate a YDMF YAML file, returning a plain Python object.

    Raises ``strictyaml.YAMLValidationError`` on schema violations.
    """
    text = Path(path).read_text()
    parsed = strictyaml.load(text, YDMF_SCHEMA)
    return parsed.data


def validate_ydmf(text: str) -> AnyType:
    """Validate a YDMF document already loaded as a string; returns plain data."""
    parsed = strictyaml.load(text, YDMF_SCHEMA)
    return parsed.data
