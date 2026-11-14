"""Placement ▸ PlanTUS: choose the transducer position on PlanTUS's maps.

PlanTUS (Lueckel et al., 2025) scores every skin vertex on five criteria. A run computes
them (the Automatic button: headless, places the composite-score optimum) and writes
plantus_metrics.mat to the placement run folder. This page draws the skin coloured by a
criterion or the composite score; a click picks a vertex, shows its values and the
transducer there, and Place here sets placement.plantus.vertex and places it without
recomputing the maps. The run itself goes through the Run tab (MATLAB), as the other
PlanTUS actions do.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import QObject, Qt, QUrl, pyqtSignal, pyqtSlot
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QHBoxLayout, QLabel, QMessageBox, QPushButton,
                             QVBoxLayout, QWidget)

from .. import mesh3d as m3
from .. import plantus as pt
from . import style
from .empty_state import EmptyStack
from .view3d import HAVE_Q3D, can_render, from_scene, look_rotation, to_scene

if HAVE_Q3D:
    from PyQt6.QtQuickWidgets import QQuickWidget

    from .view3d import MeshGeometry

TASK = "PlanTUS placement"
TARGET = "#1e88e5"
BEST = "#2e7d32"
PLACED = "ink"              # palette role: the placed position in the theme's ink (dark mode too)
CLICKED = "#ef6c00"
HOUSING = "#6666ff"
VIEW_DISTANCE_MM = 330.0

QML = """
import QtQuick
import QtQuick3D
import QtQuick3D.Helpers

Item {
    id: root
    property color background: "#ffffff"
    property var view: null          // [cx, cy, cz, qw, qx, qy, qz, distance]: reset the camera
    onViewChanged: {
        if (!view) return
        origin.position = Qt.vector3d(view[0], view[1], view[2])
        origin.rotation = Qt.quaternion(view[3], view[4], view[5], view[6])
        camera.position = Qt.vector3d(0, 0, view[7])
        camera.eulerRotation = Qt.vector3d(0, 0, 0)
    }
    View3D {
        id: view3d
        anchors.fill: parent
        environment: SceneEnvironment {
            clearColor: root.background
            backgroundMode: SceneEnvironment.Color
            antialiasingMode: SceneEnvironment.MSAA
            antialiasingQuality: SceneEnvironment.High
        }
        Node {
            id: origin
            PerspectiveCamera {
                id: camera
                position: Qt.vector3d(0, 0, 380)
                clipNear: 1
                clipFar: 5000
                DirectionalLight { brightness: 0.9 }
            }
        }
        DirectionalLight { eulerRotation.x: -60; brightness: 0.35 }
        Repeater3D {
            model: skin
            Model {
                geometry: modelData
                materials: DefaultMaterial {
                    vertexColorsEnabled: true
                    diffuseColor: "white"
                    cullMode: Material.NoCulling
                }
            }
        }
        Repeater3D {
            model: meshes
            Model {
                geometry: modelData
                materials: PrincipledMaterial {
                    baseColor: modelData.color
                    opacity: modelData.opacity
                    roughness: 0.65
                    metalness: 0.0
                    cullMode: Material.NoCulling
                    alphaMode: modelData.opacity < 1 ? PrincipledMaterial.Blend : PrincipledMaterial.Opaque
                    depthDrawMode: modelData.opacity < 1 ? Material.NeverDepthDraw : Material.OpaqueOnlyDepthDraw
                }
            }
        }
    }
    OrbitCameraController {
        anchors.fill: parent; origin: origin; camera: camera
        TapHandler {
            acceptedButtons: Qt.LeftButton | Qt.RightButton
            onTapped: function(eventPoint, button) {
                var p = eventPoint.position
                var a = view3d.mapTo3DScene(Qt.vector3d(p.x, p.y, 0))
                var b = view3d.mapTo3DScene(Qt.vector3d(p.x, p.y, 1000))
                picker.pick(a.x, a.y, a.z, b.x, b.y, b.z)
            }
        }
    }
}
"""


class PickBridge(QObject):
    """QML calls pick() with two scene points on the ray under the click."""

    picked = pyqtSignal(object, object)      # origin, direction (RAS mm)

    @pyqtSlot(float, float, float, float, float, float)
    def pick(self, ax, ay, az, bx, by, bz) -> None:
        a, b = from_scene([[ax, ay, az], [bx, by, bz]])
        if np.linalg.norm(b - a) > 0:
            self.picked.emit(a, b - a)


def run_placement(window, mode: str, task: str, interactive: bool = False) -> str | None:
    """Run only the placement stage (PlanTUS, TUS_entry, heuristic): a job where the Run tab
    submits to a cluster, else in MATLAB here; PlanTUS's Workbench picking needs a display and
    always runs here. None if nothing started (partial_run.run_partial)."""
    from .. import launcher
    from .partial_run import run_partial
    if window.model.get("placement.mode") != mode:
        QMessageBox.information(window, task, f"Set placement.mode to '{mode}' first.")
        return None
    return run_partial(window, task, mode, f"the {mode} placement",
                       launcher.placement_overrides(mode, interactive), launcher.PLACEMENT_ONLY,
                       local_only=interactive)


def run_plantus(window, interactive: bool, vertex: int | None) -> str | None:
    """PlanTUS placement only: vertex None computes the maps (and places the optimum, or the
    vertex picked in Workbench); a vertex places it on the maps of the earlier run.
    placement.plantus.vertex is saved with the config, so the study keeps the choice."""
    if window.model.get("placement.mode") == "plantus":
        window.model.set("placement.plantus.vertex", [] if vertex is None else int(vertex))
    return run_placement(window, "plantus", TASK, interactive)


class PlantusPage(QWidget):
    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self.maps: pt.PlantusMaps | None = None
        self._loaded_key = None                 # (path, mtime) of the maps shown
        self.tx: dict = {}
        self.clicked: int | None = None
        self._quick = None
        self._bridge = PickBridge()
        self._bridge.picked.connect(self._on_pick)
        self._skin_items: list = []
        self._items: list = []
        self._skin_key = None
        self._skin_verts: np.ndarray | None = None   # the skin as drawn (face scrambled)
        self._view_for = None
        self._face = (None, None)               # (segmentation, face region)
        self._waiting = False
        self._watching = False
        self._stale = False

        v = QVBoxLayout(self)
        intro = QLabel("PlanTUS scores every skin position for this target. Colour the skin by the "
                       "composite score or one criterion, click a position to see its values and the "
                       "transducer there, then <b>Place here</b>. Drag to turn, scroll to zoom.")
        intro.setObjectName("CardSub")
        intro.setWordWrap(True)
        v.addWidget(intro)

        bar = QHBoxLayout()
        bar.addWidget(QLabel("Map"))
        self.map_box = QComboBox()
        for key, label, _u, _b in pt.MAPS:
            self.map_box.addItem(label, key)
        self.map_box.currentIndexChanged.connect(lambda _i: self._rebuild(skin=True))
        bar.addWidget(self.map_box)
        self.show_face = QCheckBox("Face", checked=False)
        self.show_face.setToolTip("Draw the face as segmented. Off by default: the skin in front of and "
                                  "below the eyes is scrambled, so the view does not identify the participant "
                                  "(the face is excluded from placement anyway).")
        self.show_face.toggled.connect(lambda _on: self._rebuild(skin=True))
        bar.addWidget(self.show_face)
        bar.addStretch(1)
        self.compute_btn = QPushButton("Compute maps")
        self.compute_btn.setToolTip("Run PlanTUS for this subject and target: computes the maps and places "
                                    "the composite-score optimum (clears placement.plantus.vertex)")
        self.compute_btn.clicked.connect(self.compute)
        bar.addWidget(self.compute_btn)
        reset = QPushButton("Reset view")
        reset.clicked.connect(lambda: self._reset_view(force=True))
        bar.addWidget(reset)
        v.addLayout(bar)

        self._content = QWidget()
        body = QHBoxLayout(self._content)
        body.setContentsMargins(0, 0, 0, 0)
        self._view_box = QVBoxLayout()
        body.addLayout(self._view_box, 1)
        side = QVBoxLayout()
        self.readout = QLabel()
        self.readout.setWordWrap(True)
        self.readout.setTextFormat(Qt.TextFormat.RichText)
        self.readout.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.readout.setMinimumWidth(270)
        side.addWidget(self.readout)
        self.place_btn = QPushButton("Place here")
        self.place_btn.setToolTip("Set placement.plantus.vertex to this position and place the transducer "
                                  "there (Run tab; the maps are not recomputed)")
        self.place_btn.clicked.connect(self.place)
        side.addWidget(self.place_btn)
        self.best_btn = QPushButton("Go to the optimum")
        self.best_btn.setToolTip("Select the position with the best composite score")
        self.best_btn.clicked.connect(lambda: self.select(self.maps.best if self.maps else None))
        side.addWidget(self.best_btn)
        side.addStretch(1)
        self.status = QLabel()
        self.status.setObjectName("CardSub")
        self.status.setWordWrap(True)
        side.addWidget(self.status)
        body.addLayout(side)
        self.stack = EmptyStack(self._content)
        v.addWidget(self.stack, 1)
        self.legend = QLabel()
        self.legend.setObjectName("CardSub")
        self.legend.setWordWrap(True)
        v.addWidget(self.legend)

    @property
    def model(self):
        return self.window.model

    # ── refresh ─────────────────────────────────────────────────────────
    def refresh(self) -> None:
        if self.window.isVisible() and not self.isVisible():
            self._stale = True
            return
        is_mode = self.model.get("placement.mode") == "plantus"
        self.compute_btn.setEnabled(is_mode)
        path = pt.metrics_path(self.model)
        if path is None:
            self.stack.show_empty("No placement run folder yet",
                                  "The PlanTUS maps live in the placement run folder, which needs path.sim "
                                  "(or path.placement) and the subject; set them on Setup ▸ Paths.")
            return
        try:
            key = (path, path.stat().st_mtime)
        except OSError:
            key = None
        if key is None:
            self.maps = None
            self._loaded_key = None
            body = ("Compute maps runs PlanTUS for this subject and target (a few minutes; as a cluster "
                    "job when the Run tab submits there) and places the composite-score optimum; the "
                    "maps then appear here." if is_mode else
                    "Set placement.mode to 'plantus' (Placement tab) and the PlanTUS settings first.")
            self.stack.show_empty("No PlanTUS maps yet", f"{body}<br>Looked for {path}",
                                  "Compute maps" if is_mode else None, self.compute if is_mode else None)
            return
        if key != self._loaded_key:
            try:
                self.maps = pt.load(path)
            except Exception as exc:          # noqa: BLE001 — shown, not raised
                self.maps = None
                self._loaded_key = None
                self.stack.show_empty("Could not read the PlanTUS maps", f"{path}: {exc}")
                return
            self._loaded_key = key
            self.tx = pt.read_tx_config(path.parent)
            self._skin_key = None
            chosen = self.model.get("placement.plantus.vertex")
            self.clicked = self.maps.selected if self.maps.selected is not None else self.maps.best
            if isinstance(chosen, (int, float)) and 0 <= int(chosen) < self.maps.n:
                self.clicked = int(chosen)
        if not can_render():
            self.stack.show_empty("The 3D view needs a display",
                                  "Qt Quick 3D could not be used here (no GPU display, or PyQt6's "
                                  "QtQuick3D module is missing).")
            return
        self.stack.show_content()
        self._ensure_quick()
        self._rebuild(skin=True)
        self._reset_view()

    def showEvent(self, event) -> None:        # noqa: N802
        super().showEvent(event)
        if self._stale:
            self._stale = False
            self.refresh()

    def _ensure_quick(self) -> None:
        if self._quick is not None:
            return
        from ..localcopy import cache_root
        qml = cache_root() / "plantus_3d.qml"            # local disk: QML must load from a file
        qml.parent.mkdir(parents=True, exist_ok=True)
        if not qml.is_file() or qml.read_text(encoding="utf-8") != QML:
            qml.write_text(QML, encoding="utf-8")
        self._quick = QQuickWidget()
        self._quick.setResizeMode(QQuickWidget.ResizeMode.SizeRootObjectToView)
        self._quick.setMinimumHeight(320)
        ctx = self._quick.rootContext()
        ctx.setContextProperty("skin", [])
        ctx.setContextProperty("meshes", [])
        ctx.setContextProperty("picker", self._bridge)
        self._quick.setSource(QUrl.fromLocalFile(str(qml)))
        self._view_box.addWidget(self._quick, 1)
        from .export import attach_grab, subject_stem
        self.export_button = attach_grab(self._quick, "plantus_3d",
                                         subject_stem(lambda: self.model, "plantus_3d"))

    # ── scene ───────────────────────────────────────────────────────────
    def map_key(self) -> str:
        return self.map_box.currentData() or "score"

    def skin_vertices(self) -> np.ndarray:
        """The skin as drawn: the face scrambled unless Face is ticked (needs the segmentation
        the Viewer loads; without one the skin is drawn as is)."""
        verts = self.maps.coords
        if self.show_face.isChecked():
            return verts
        region = self._face_region()
        return verts if region is None else m3.scramble_region(verts, self.maps.faces, *region)[0]

    def _face_region(self):
        from .main_window import VIEWER
        viewer = self.window.pages.get(VIEWER)
        seg = None
        if viewer is not None:
            try:
                viewer._load_images()
                seg = viewer.seg
            except Exception:                 # noqa: BLE001 — no segmentation: draw unscrambled
                seg = None
        if seg is None:
            return None
        if self._face[0] is not seg:
            step = max(1, int(round(2.0 / float(np.min(seg.voxel_size)))))
            sub = seg.data[::step, ::step, ::step]
            self._face = (seg, m3.face_region(sub, seg.affine @ np.diag([step, step, step, 1.0])))
        return self._face[1]

    def marker_meshes(self) -> list[tuple]:
        """(vertices, faces, colour, opacity): target, optimum, placed vertex, the clicked
        vertex with the transducer PlanTUS would put there."""
        m = self.maps
        out = []
        if np.all(np.isfinite(m.target)):
            out.append((*m3.sphere(m.target, 2.5), TARGET, 1.0))
        if m.best is not None:
            out.append((*m3.sphere(m.coords[m.best], 2.2), BEST, 1.0))
        if m.selected is not None and m.selected != m.best:
            out.append((*m3.sphere(m.coords[m.selected], 2.2), style.active()[PLACED], 1.0))
        if self.clicked is not None:
            p = pt.preview(m, self.clicked, self.tx)
            out.append((*m3.sphere(p["skin"], 2.0), CLICKED, 1.0))
            out.append((*m3.cylinder(p["centre"], m.target, 0.6), CLICKED, 1.0))
            if p["diameter"] > 0:
                plane = float(self.tx.get("plane_offset", 0.0) or 0.0)
                d, depth, front, _known = m3.housing_size(p["diameter"], plane, self.housing())
                out.append((*m3.housing_puck(d, depth, front, p["centre"], p["axis"]), HOUSING, 0.55))
        return out

    def housing(self) -> dict:
        """The housing block of the transducer PlanTUS places (transducer(1)), as Placement ▸ 3D
        draws it; {} without one."""
        model = self.model
        lib = getattr(model, "library", None)
        try:
            txs = model.effective_transducers()
        except Exception:                     # noqa: BLE001 — no housing rather than no preview
            return {}
        return lib.housing(txs[0].get("serial")) if lib is not None and txs else {}

    def _rebuild(self, skin: bool = False) -> None:
        if self._quick is None or self.maps is None:
            return
        ctx = self._quick.rootContext()
        key = (self.map_key(), self.show_face.isChecked(), self._loaded_key)
        if skin and key != self._skin_key:
            rgba = pt.colours(pt.goodness(self.maps, self.map_key()))
            if self._skin_key is None or self._skin_key[1:] != key[1:]:
                self._skin_verts = self.skin_vertices()
            self._skin_items = [MeshGeometry(self._skin_verts, self.maps.faces, "#ffffff", 1.0, rgba)]
            ctx.setContextProperty("skin", self._skin_items)
            self._skin_key = key
        meshes = sorted(self.marker_meshes(), key=lambda mm: mm[3] < 1)
        self._items = [MeshGeometry(v, f, c, a) for v, f, c, a in meshes]
        ctx.setContextProperty("meshes", self._items)
        root = self._quick.rootObject()
        if root is not None:
            root.setProperty("background", style.active()["surface"])
        self._update_side()

    def _update_side(self) -> None:
        m = self.maps
        if self.clicked is None:
            self.readout.setText("Click the skin to choose a position.")
        else:
            rows = "".join(f"<tr><td>{k}</td><td style='padding-left:8px'><b>{v}</b></td></tr>"
                           for k, v in m.readout(self.clicked))
            self.readout.setText(f"<table>{rows}</table>")
        is_mode = self.model.get("placement.mode") == "plantus"
        self.place_btn.setEnabled(is_mode and self.clicked is not None and not self._waiting)
        self.best_btn.setEnabled(m.best is not None)
        chosen = self.model.get("placement.plantus.vertex")
        placed = ("placed: vertex " + str(m.selected) + f" ({m.selection_mode})") if m.selected is not None else ""
        cfg = (f"placement.plantus.vertex = {int(chosen)}" if isinstance(chosen, (int, float))
               else "placement.plantus.vertex = [] (optimum)")
        mode_note = "" if is_mode else "<br>placement.mode is not 'plantus': placing is off."
        self.status.setText(f"{m.path.parent}<br>{placed}<br>{cfg}{mode_note}")
        key, label = self.map_key(), pt.MAP_LABELS[self.map_key()]
        vals = m.maps[key][np.isfinite(m.maps[key]) & (m.avoidance != 0)]
        unit = {k: u for k, _l, u, _b in pt.MAPS}[key]
        rng = f" ({vals.min():.1f}–{vals.max():.1f} {unit})".replace(" )", ")") if len(vals) else ""
        better = ("closer to the optimal distance " f"{m.optimal_distance:.0f} mm"
                  if key == "skin_target_distances" else
                  {"low": "lower", "high": "higher"}[{k: b for k, _l, _u, b in pt.MAPS}[key]])
        weights = ", ".join(f"{k.replace('_', ' ')} {w:g}" for k, w in m.weights.items())
        self.legend.setText(
            f"{label}{rng}: bright = better ({better}); grey = excluded (avoidance area or no value) · "
            f"<span style='color:{TARGET}'>●</span> target · <span style='color:{BEST}'>●</span> optimum · "
            f"<span style='color:{style.active()[PLACED]}'>●</span> placed · <span style='color:{CLICKED}'>●</span> chosen, "
            f"with the beam axis and the transducer · click to choose"
            + (f" · weights: {weights}" if weights else ""))

    def _reset_view(self, force: bool = False) -> None:
        """Centre on the target, seen from the chosen (else optimum) position."""
        if self._quick is None or self._quick.rootObject() is None or self.maps is None:
            return
        key = self._loaded_key
        if not force and self._view_for == key:
            return
        self._view_for = key
        m = self.maps
        centre = m.target if np.all(np.isfinite(m.target)) else m.coords.mean(axis=0)
        v = self.clicked if self.clicked is not None else m.best
        outward = (m.coords[v] - centre) if v is not None else np.array([1.0, 0.6, 0.4])
        if np.linalg.norm(outward) == 0:
            outward = np.array([1.0, 0.6, 0.4])
        q = look_rotation(to_scene(outward / np.linalg.norm(outward))[0])
        c = to_scene(centre)[0]
        self._quick.rootObject().setProperty(
            "view", [float(c[0]), float(c[1]), float(c[2]), *q, VIEW_DISTANCE_MM])

    # ── interaction ─────────────────────────────────────────────────────
    def _on_pick(self, origin, direction) -> None:
        if self.maps is None:
            return
        verts = self._skin_verts if self._skin_verts is not None else self.maps.coords
        v = pt.pick_vertex(origin, direction, verts, self.maps.faces)
        if v is not None:
            self.select(v)

    def select(self, v: int | None) -> None:
        if v is None or self.maps is None:
            return
        self.clicked = int(v)
        self._rebuild()

    def place(self) -> None:
        if self.maps is None or self.clicked is None:
            return
        if not self.maps.eligible(self.clicked):
            ok = QMessageBox.question(
                self, TASK, "This position is not eligible under the composite score (see the readout). "
                "Place the transducer there anyway?")
            if ok != QMessageBox.StandardButton.Yes:
                return
        self._start(False, self.clicked)

    def compute(self) -> None:
        self._start(False, None)

    def _start(self, interactive: bool, vertex: int | None) -> None:
        from .main_window import RUN
        how = run_plantus(self.window, interactive, vertex)
        run = self.window.pages[RUN]
        if how == "local" and not self._waiting:
            self._waiting = True
            run.finished.connect(self._done)
        elif how == "job" and not self._watching:
            self._watching = True                   # reload the maps once the job has written them
            run.jobs_changed.connect(self._jobs_changed)
        if how is not None and self.maps is not None:
            self._update_side()

    def _jobs_changed(self) -> None:
        if self.isVisible():
            self.refresh()                          # a stat; reloads only when the maps changed

    def _done(self, task: str, ok: bool) -> None:
        from .main_window import RUN
        if task != TASK:
            return
        self.window.pages[RUN].finished.disconnect(self._done)
        self._waiting = False
        self.window.statusBar().showMessage(
            "PlanTUS finished: the Viewer shows the placement." if ok else
            "PlanTUS failed: see the Run tab's log.", 10000)
        self.refresh()
