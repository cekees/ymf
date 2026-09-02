import pytest

from ymf.data_sources import (
    DEFAULT_MANNING_N_TABLE,
    ManningLookupError,
    KNOWN_DATA_SOURCES,
    load_manning_n_table,
    manning_n_for_land_cover,
    resolve_roughness_source,
    validate_source_record,
)


def test_manning_n_for_land_cover_known():
    n = manning_n_for_land_cover("forest_dense")
    assert n == DEFAULT_MANNING_N_TABLE["forest_dense"]
    assert 0.10 <= n <= 0.25  # per report §4.1 range for dense forest


def test_manning_n_for_land_cover_unknown_raises():
    with pytest.raises(ManningLookupError):
        manning_n_for_land_cover("not_a_real_land_cover_class")


def test_manning_n_for_nlcd_class():
    # NLCD 42 = evergreen forest -> forest_dense
    n = manning_n_for_land_cover("unused", nlcd_class=42)
    assert n == DEFAULT_MANNING_N_TABLE["forest_dense"]


def test_manning_n_for_unmapped_nlcd_class_raises():
    with pytest.raises(ManningLookupError):
        manning_n_for_land_cover("unused", nlcd_class=999)


def test_load_manning_n_table_default():
    table = load_manning_n_table()
    assert table == DEFAULT_MANNING_N_TABLE


def test_load_manning_n_table_csv_override(tmp_path):
    csv_path = tmp_path / "manning_n_lookup.csv"
    csv_path.write_text("land_cover,n\nforest_dense,0.20\ncustom_class,0.05\n")
    table = load_manning_n_table(csv_path)
    assert table["forest_dense"] == 0.20  # overridden
    assert table["custom_class"] == 0.05  # added
    assert table["open_water"] == DEFAULT_MANNING_N_TABLE["open_water"]  # untouched


def test_resolve_roughness_source_land_cover_key():
    roughness_source = {
        "type": "land_cover_mapping",
        "land_cover_data": {"source": "ESA CCI", "resolution": "300 m"},
        "reference": "Arcement and Schneider (1989)",
    }
    n = resolve_roughness_source(roughness_source, "cropland")
    assert n == DEFAULT_MANNING_N_TABLE["cropland"]


def test_resolve_roughness_source_nlcd_int_class():
    roughness_source = {
        "type": "land_cover_mapping",
        "land_cover_data": {"source": "USGS NLCD", "resolution": "30 m"},
    }
    n = resolve_roughness_source(roughness_source, 41)  # deciduous forest
    assert n == DEFAULT_MANNING_N_TABLE["forest_dense"]


def test_validate_source_record_known_source_no_warnings():
    record = {"source": "CHIRPS", "source_url": KNOWN_DATA_SOURCES["CHIRPS"].url}
    assert validate_source_record(record) == []


def test_validate_source_record_unknown_source_warns():
    record = {"source": "SomeRandomDataset", "source_url": "https://example.com"}
    warnings = validate_source_record(record)
    assert len(warnings) == 1
    assert "not in the known data source catalog" in warnings[0]


def test_validate_source_record_mismatched_url_warns():
    record = {"source": "CHIRPS", "source_url": "https://mirror.example.com/chirps"}
    warnings = validate_source_record(record)
    assert len(warnings) == 1
    assert "does not match the catalog URL" in warnings[0]


def test_validate_source_record_no_source_is_fine():
    assert validate_source_record({}) == []
