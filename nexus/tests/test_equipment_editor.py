"""Equipment editor rules (prestus_nexus.equipment_editor), as calibration_gui's 'New…' form."""
from pathlib import Path
from types import SimpleNamespace

from prestus_nexus import equipment_editor as ed
from prestus_nexus import yamlio


def _form(**kw):
    f = ed.TransducerForm(serial="IS_TEST_001", name="", manufact="Imasonic", freq_khz=300, curv_radius_mm=75,
                       dist_geom_ep_mm=64, min_foc=30, max_foc=80, n_channels=2,
                       elem_id_mm="0, 30, 40, 50", elem_od_mm="29, 39, 49, 60")
    for k, v in kw.items():
        setattr(f, k, v)
    return f


def test_parse_list():
    assert ed.parse_list("0, 30.2 [42.1; 51]") == [0.0, 30.2, 42.1, 51.0]
    assert ed.parse_list("1, x") is None
    assert ed.parse_list("") == []


def test_valid_entry(tmp_path):
    entry, issues = ed.build_entry(_form(), tmp_path)
    assert issues == []
    assert entry["type"] == "transducer" and entry["serial"] == "IS_TEST_001"
    assert entry["name"] == "IS_TEST_001"                  # name defaults to the serial
    assert entry["n_elem"] == 2
    a = entry["transducer"]["annular"]
    assert a["elem_n"] == 4 and a["elem_od_mm"][-1] == 60.0
    assert entry["transducer"]["freq_hz"] == 300e3
    assert "combos" not in entry


def test_issues():
    _, issues = ed.build_entry(_form(serial="1bad", min_foc=90, n_channels=3, elem_od_mm="29, 39, 39, 60"),
                               Path("."))
    text = " ".join(issues)
    assert "Serial must start with a letter" in text
    assert "minimum must be below maximum" in text
    assert "outer diameter must exceed" in text
    assert "multiple of the physical channel count" in text
    _, issues = ed.build_entry(_form(elem_id_mm="0, 30", elem_od_mm="29, 39, 49"), Path("."))
    assert any("2 inner but 3 outer" in i for i in issues)


def test_combos_profile_path_and_new_ds(tmp_path):
    root = tmp_path / "local" / "profiles"
    csv = root / "IS_TEST_001" / "tran-IS_TEST_001_ds-NEWDS_axial.csv"
    csv.parent.mkdir(parents=True)
    csv.write_text("x\n")
    new = ed.NewDrivingSystem("NEWDS", "New driver", "Acme", 8)
    form = _form(combos=[ed.Combo("NEWDS", str(csv), new), ed.Combo("SC_203_035")])
    entry, issues = ed.build_entry(form, root)
    assert issues == []
    assert entry["combos"]["NEWDS"] == {"char_data_path": "IS_TEST_001/tran-IS_TEST_001_ds-NEWDS_axial.csv"}
    assert entry["combos"]["SC_203_035"] == {}
    files = ed.target_files(entry, tmp_path / "lib", None, [new])
    names = sorted(p.name for p in files)
    assert names == ["ds-NEWDS_driver.yaml", "tran-IS_TEST_001_transducer.yaml"]
    ed.write(files)
    ds = yamlio.load_file(tmp_path / "lib" / "ds-NEWDS_driver.yaml")
    assert ds == {"type": "driving_system", "serial": "NEWDS", "name": "New driver", "manufact": "Acme",
                  "available_ch": 8}
    back = yamlio.load_file(tmp_path / "lib" / "tran-IS_TEST_001_transducer.yaml")
    assert back["serial"] == "IS_TEST_001" and back["combos"]["NEWDS"]["char_data_path"].endswith("_axial.csv")


def test_combo_rules():
    _, issues = ed.build_entry(_form(combos=[ed.Combo("", "a.csv")]), Path("."))
    assert any("belongs to a driving system" in i for i in issues)
    _, issues = ed.build_entry(_form(combos=[ed.Combo("A"), ed.Combo("A")]), Path("."))
    assert any("listed twice" in i for i in issues)
    _, issues = ed.build_entry(_form(n_channels=4, combos=[ed.Combo("D", "", ed.NewDrivingSystem("D", available_ch=2))]),
                               Path("."))
    assert any("fewer than the transducer's 4" in i for i in issues)


def test_edit_existing_keeps_fields_and_combos(tmp_path):
    raw = {"type": "transducer", "serial": "CTX_X", "name": "CTX X", "manufact": "SC", "n_elem": 4,
           "min_foc": 10, "max_foc": 60, "housing": {"diameter_mm": 70},
           "transducer": {"annular": {"elem_n": 4, "elem_id_mm": [0, 30, 40, 50], "elem_od_mm": [29, 39, 49, 60],
                                      "curv_radius_mm": 63, "dist_geom_ep_mm": 52}, "freq_hz": 250e3},
           "combos": {"SC_203_035": {"char_data_path": "CTX_X/p.csv", "drive_limit": {"max_elem_p_amp_pa": 1}}}}
    dev = SimpleNamespace(serial="CTX_X", name="CTX X", manufact="SC", model="", freq_hz=250e3, min_foc=10,
                          max_foc=60, n_elem=4, geometry=raw["transducer"], user=True,
                          file=str(tmp_path / "own" / "tran-CTX_X_transducer.yaml"))
    form, base = ed.form_from_device(dev, raw, keep_serial=True)
    assert form.serial == "CTX_X" and [c.ds_serial for c in form.combos] == ["SC_203_035"]
    form.combos.append(ed.Combo("IGT_32", ""))
    entry, issues = ed.build_entry(form, tmp_path, base)
    assert issues == []
    assert entry["housing"] == {"diameter_mm": 70}                          # kept from the raw entry
    assert entry["combos"]["SC_203_035"]["drive_limit"] == {"max_elem_p_amp_pa": 1}
    assert entry["combos"]["SC_203_035"]["char_data_path"] == "CTX_X/p.csv"
    assert "IGT_32" in entry["combos"]
    files = ed.target_files(entry, tmp_path / "elsewhere", dev, [])
    assert list(files) == [Path(dev.file)]                                  # rewrites its own file
    copy_form, copy_base = ed.form_from_device(dev, raw, keep_serial=False)
    assert copy_form.serial == "" and copy_form.combos == [] and copy_base is None


def test_matrix_entry(tmp_path):
    form = _form(kind="matrix", curv_radius_mm=62.94, dist_geom_ep_mm=52.38, elem_id_mm="", elem_od_mm="")
    form.matrix = ed.MatrixFields(grid="fibonacci", fib_n=256, fib_kerf_mm=0.5, outer_diameter_mm=60)
    entry, issues = ed.build_entry(form, tmp_path)
    assert any("layout not computed" in i for i in issues)            # count comes from MATLAB
    entry, issues = ed.build_entry(form, tmp_path, matrix_n_elem=256)
    assert issues == []
    m = entry["transducer"]["matrix"]
    assert entry["n_elem"] == 256 and "annular" not in entry["transducer"]
    assert m["matrix_shape"]["define_here"]["grid_shape"] == {"type": "fibonacci",
                                                             "fibonacci": {"elem_n": 256, "kerf_mm": 0.5}}
    form.matrix.outer_diameter_mm = 130                               # radius 62.94 < 65
    _, issues = ed.build_entry(form, tmp_path, matrix_n_elem=256)
    assert any("half the aperture" in i for i in issues)
    _, issues = ed.build_entry(form, tmp_path, matrix_n_elem=0)
    assert any("No element" in i for i in issues)


def test_matrix_round_trip_from_device():
    geo = {"matrix": {"steering": "1D", "elem_shape": "disc", "elem_width_mm": 2.0, "elem_height_mm": 2.0,
                      "outer_diameter_mm": 50.0, "is_curved": False, "curv_radius_mm": 70.0,
                      "dist_geom_ep_mm": 60.0,
                      "matrix_shape": {"type": "define_here", "define_here": {"grid_shape": {
                          "type": "rect", "rect": {"elem_n_row": 8, "elem_n_col": 4, "elem_spacing_height_mm": 1.0,
                                                   "elem_spacing_width_mm": 0.5, "sparsity_factor": 1}}}}},
           "freq_hz": 500e3}
    dev = SimpleNamespace(serial="MX", name="MX", manufact="", model="", freq_hz=500e3, min_foc=20, max_foc=80,
                          n_elem=32, geometry=geo, user=False, file="")
    form, _ = ed.form_from_device(dev, None, keep_serial=False)
    assert form.kind == "matrix" and form.matrix.grid == "rect" and (form.matrix.n_row, form.matrix.n_col) == (8, 4)
    form.serial = "MX_COPY"
    entry, issues = ed.build_entry(form, Path("."), matrix_n_elem=32)
    assert issues == [] and entry["transducer"]["matrix"]["matrix_shape"] == geo["matrix"]["matrix_shape"]


def test_shipped_targets(tmp_path):
    shipped = tmp_path / "equipment"
    (shipped / "model-CTX_500").mkdir(parents=True)
    entry = {"serial": "CTX_500_099", "model": "CTX_500"}
    files = ed.target_files(entry, tmp_path / "own", None, [ed.NewDrivingSystem("D1")], shipped=shipped)
    assert sorted(p.relative_to(shipped).as_posix() for p in files) == ["ds/ds-D1_driver.yaml",
                                                                    "model-CTX_500/tran-CTX_500_099_transducer.yaml"]
    files = ed.target_files({"serial": "NEW_X"}, tmp_path / "own", None, [], shipped=shipped)
    assert list(files) == [shipped / "NEW_X.yaml"]
