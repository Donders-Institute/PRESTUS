"""Placement ▸ 3D: the head, the brain and every transducer's geometry in 3D.

The bowl is drawn as it is built in the simulation (one curved annulus per element ring, flat
for a planar device), with the exit plane, the beam axis, the target, the geometric focus and
the point the focal depth setting lands on. Drag to turn, scroll to zoom.

Rendered with Qt Quick 3D (part of PyQt6, no further dependency). The page reads images,
positions and geometry from the placement viewer (PlacementPage), as the Beam axis tab does.
Meshes are made in Python (mesh3d.py) and handed to QML as a list of geometry objects.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import QByteArray, QObject, Qt, QUrl, pyqtSlot
from PyQt6.QtGui import QGuiApplication, QVector3D
from PyQt6.QtWidgets import (QCheckBox, QHBoxLayout, QLabel, QPushButton, QSlider, QVBoxLayout,
                             QWidget)

from .. import mesh3d as m3
from .. import placement as pl
from . import style
from .empty_state import EmptyStack

try:
    from PyQt6.QtCore import pyqtProperty
    from PyQt6.QtQuick3D import QQuick3DGeometry
    from PyQt6.QtQuickWidgets import QQuickWidget
    HAVE_Q3D = True
except ImportError:          # pragma: no cover
    HAVE_Q3D = False

SCALP = ("#e8b49a", 0.30)                    # colour, default opacity (slider)
BRAIN = ("#cfc8dc", 0.60)                    # see-through by default: the target lies inside
MASK = ("#43a047", 0.55)                     # target mask of a PlanTUS / TUS_entry run
RINGS_SELECTED = ("#c62828", "#ef6c00")      # alternate rings, so each element is visible
RINGS_OTHER = ("#9e9e9e", "#bdbdbd")
EXIT_PLANE = ("#ffd54f", 0.45)
HOUSING = "#6666ff"                          # BabelBrain's transducer colour (0.4, 0.4, 1.0)
HOUSING_OTHER = "#b0bec5"
TARGET = "#1e88e5"
GEOM_FOCUS = "#757575"
SETTING_FOCUS = "#ff9800"
VIEW_DISTANCE_MM = 260.0

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
                DirectionalLight { brightness: 0.9 }     // a headlight: lit from where you look
            }
        }
        DirectionalLight { eulerRotation.x: -60; brightness: 0.35 }
        Repeater3D {
            model: meshes
            Model {
                geometry: modelData
                property string tip: modelData.tip
                pickable: tip !== ""            // only the transducers: picking stays cheap
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
    OrbitCameraController { anchors.fill: parent; origin: origin; camera: camera }
    // Hover a transducer: its name. Picked once the pointer rests, not on every move.
    HoverHandler {
        id: hover
        onPointChanged: tipTimer.restart()
        onHoveredChanged: if (!hovered) tipSink.show("")
    }
    Timer {
        id: tipTimer
        interval: 120
        onTriggered: {
            if (!hover.hovered) return
            const hit = view3d.pick(hover.point.position.x, hover.point.position.y).objectHit
            tipSink.show(hit && hit.tip ? hit.tip : "")
        }
    }
}
"""


def to_scene(p: np.ndarray) -> np.ndarray:
    """RAS mm -> Qt Quick 3D scene (y up): x = R, y = S, z = -A (a rotation, no mirror)."""
    p = np.atleast_2d(np.asarray(p, float))
    return np.c_[p[:, 0], p[:, 2], -p[:, 1]]


def from_scene(p: np.ndarray) -> np.ndarray:
    """Inverse of to_scene: Qt Quick 3D scene -> RAS mm."""
    p = np.atleast_2d(np.asarray(p, float))
    return np.c_[p[:, 0], -p[:, 2], p[:, 1]]


def look_rotation(direction: np.ndarray) -> tuple[float, float, float, float]:
    """Quaternion (w, x, y, z) turning the scene's +z onto ``direction`` (scene coordinates)
    with superior kept up, so a camera at +z behind the origin looks along -direction."""
    z = np.asarray(direction, float)
    z = z / np.linalg.norm(z)
    up = np.array([0.0, 1.0, 0.0]) if abs(z[1]) < 0.95 else np.array([0.0, 0.0, -1.0])
    x = np.cross(up, z)
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    m = np.stack([x, y, z], axis=1)
    w = np.sqrt(max(0.0, 1 + m[0, 0] + m[1, 1] + m[2, 2])) / 2
    if w > 1e-6:
        return (w, (m[2, 1] - m[1, 2]) / (4 * w), (m[0, 2] - m[2, 0]) / (4 * w),
                (m[1, 0] - m[0, 1]) / (4 * w))
    i = int(np.argmax(np.diag(m)))                  # 180° turn: largest diagonal element
    j, k = (i + 1) % 3, (i + 2) % 3
    s = np.sqrt(max(0.0, 1 + m[i, i] - m[j, j] - m[k, k])) * 2
    q = np.zeros(4)
    q[0] = (m[k, j] - m[j, k]) / s
    q[1 + i] = s / 4
    q[1 + j] = (m[j, i] + m[i, j]) / s
    q[1 + k] = (m[k, i] + m[i, k]) / s
    return tuple(float(c) for c in q)


if HAVE_Q3D:
    class MeshGeometry(QQuick3DGeometry):
        """One triangle mesh (RAS mm) with its colour and opacity, for the QML scene;
        ``rgba`` (N, 4) colours it per vertex instead (drawn with vertex colours)."""

        def __init__(self, verts: np.ndarray, faces: np.ndarray, color: str, opacity: float = 1.0,
                     rgba: np.ndarray | None = None, tip: str = ""):
            super().__init__()
            self._color, self._opacity, self._tip = color, float(opacity), tip
            pos = to_scene(verts).astype(np.float32)
            nrm = to_scene(m3.vertex_normals(verts, faces)).astype(np.float32)
            cols = [pos, nrm] + ([np.asarray(rgba, np.float32)] if rgba is not None else [])
            self.setVertexData(QByteArray(np.hstack(cols).astype(np.float32).tobytes()))
            self.setIndexData(QByteArray(np.asarray(faces, np.uint32).tobytes()))
            self.setStride(40 if rgba is not None else 24)
            self.setPrimitiveType(QQuick3DGeometry.PrimitiveType.Triangles)
            sem, ctype = QQuick3DGeometry.Attribute.Semantic, QQuick3DGeometry.Attribute.ComponentType
            self.addAttribute(sem.PositionSemantic, 0, ctype.F32Type)
            self.addAttribute(sem.NormalSemantic, 12, ctype.F32Type)
            if rgba is not None:
                self.addAttribute(sem.ColorSemantic, 24, ctype.F32Type)
            self.addAttribute(sem.IndexSemantic, 0, ctype.U32Type)
            lo, hi = pos.min(axis=0), pos.max(axis=0)
            self.setBounds(QVector3D(*map(float, lo)), QVector3D(*map(float, hi)))

        @pyqtProperty(str, constant=True)
        def color(self) -> str:
            return self._color

        @pyqtProperty(float, constant=True)
        def opacity(self) -> float:
            return self._opacity

        @pyqtProperty(str, constant=True)
        def tip(self) -> str:
            """Hover text ('' = not pickable): the transducer a ring or housing belongs to."""
            return self._tip


class _Mesh(tuple):
    """(vertices, faces, colour, opacity) with the hover text of what it draws."""
    tip = ""


def _tipped(mesh: tuple, tip: str) -> _Mesh:
    out = _Mesh(mesh)
    out.tip = tip
    return out


class _TipSink(QObject):
    """Shows the hover text the QML scene picks (as a tooltip over the 3D widget)."""

    def __init__(self, widget):
        super().__init__(widget)
        self.widget, self._shown = widget, ""

    @pyqtSlot(str)
    def show(self, text: str) -> None:
        from PyQt6.QtGui import QCursor
        from PyQt6.QtWidgets import QApplication, QToolTip
        if QApplication.mouseButtons() != Qt.MouseButton.NoButton:
            text = ""                           # turning the view: no tooltip
        if text == self._shown:
            return
        self._shown = text
        if text:
            QToolTip.showText(QCursor.pos(), text, self.widget)
        else:
            QToolTip.hideText()


def target_mask_path(model) -> "Path | None":
    """The target mask of the placement run, if the mode has one: placement.tus_entry.target_mask,
    else the mask the run wrote to its placement run folder (TUS_entry
    sub-XXX_<target>_target_mask.nii.gz, PlanTUS sub-XXX_<target>_mask.nii.gz). Not checked
    on disk."""
    from pathlib import Path

    from ..model import is_empty
    from ..outputs import dir_placement_run
    mode = model.get("placement.mode")
    if mode not in ("tus_entry", "plantus"):
        return None
    if mode == "tus_entry" and not is_empty(model.get("placement.tus_entry.target_mask")):
        return Path(str(model.get("placement.tus_entry.target_mask")))
    sid = model.get("subject_id")
    run = dir_placement_run(model, mode)
    if run is None or not isinstance(sid, (int, float)) or isinstance(sid, bool):
        return None
    name = model.get(f"placement.{mode}.target_name")
    name = "target" if is_empty(name) else str(name)
    suffix = "target_mask" if mode == "tus_entry" else "mask"
    return run / f"sub-{int(sid):03d}_{name}_{suffix}.nii.gz"


def can_render() -> bool:
    """Qt Quick 3D needs a real display (GPU); offscreen (tests, headless) shows a note."""
    return HAVE_Q3D and QGuiApplication.platformName() not in ("offscreen", "minimal")


class Placement3DPage(QWidget):
    def __init__(self, viewer, parent: QWidget | None = None):
        super().__init__(parent)
        self.v = viewer                     # PlacementPage: images, positions, geometry
        self.model = viewer.model
        self._quick = None
        self._items: list = []              # keep the geometry objects alive while shown
        self._head: tuple | None = None     # (segmentation, {'scalp': mesh, 'brain': mesh})
        self._view_for: tuple | None = None
        self._stale = False
        self._housing_known = True          # the selected device's housing size is on record
        self._mask: tuple | None = None     # ((path, mtime), mesh) of the target mask drawn

        v = QVBoxLayout(self)
        bar = QHBoxLayout()
        self.tx_label = QLabel()
        self.tx_label.setObjectName("CardSub")
        bar.addWidget(self.tx_label)
        bar.addStretch(1)
        self.show_scalp = QCheckBox("Scalp", checked=True)
        self.show_face = QCheckBox("Face", checked=False)
        self.show_face.setToolTip("Draw the face as segmented. Off by default: the scalp in front of "
                                  "and below the eyes is scrambled (smoothed and jittered), so the head "
                                  "keeps its outline but the view does not identify the participant.")
        self.show_brain = QCheckBox("Brain", checked=True)
        self.show_mask = QCheckBox("Mask", checked=True)
        self.show_mask.setToolTip("The target mask of the placement run (PlanTUS, TUS_entry)")
        self.show_mask.setVisible(False)    # shown once a mask is found
        self.show_housing = QCheckBox("Housing", checked=True)
        self.show_housing.setToolTip("The transducer's housing as a cylinder holding the element rings, "
                                     "sized from the device's equipment entry (generic size if it has "
                                     "none). Its front face is drawn in the exit plane; a gel pad or "
                                     "membrane can leave a gap to the scalp in practice.")
        for box in (self.show_scalp, self.show_face, self.show_brain, self.show_mask, self.show_housing):
            box.toggled.connect(lambda _on: self._rebuild())
            bar.addWidget(box)
        self.scalp_opacity = self._opacity_slider(
            bar, self.show_scalp, SCALP[1], "Scalp opacity: lower it to see the brain and the target")
        self.brain_opacity = self._opacity_slider(
            bar, self.show_brain, BRAIN[1], "Brain opacity: lower it to see the target inside")
        self.mask_opacity = self._opacity_slider(bar, self.show_mask, MASK[1], "Mask opacity")
        self.mask_opacity.setVisible(False)
        self.housing_opacity = self._opacity_slider(
            bar, self.show_housing, 1.0, "Housing opacity: lower it to see the element rings and the "
            "exit plane inside")
        reset = QPushButton("Reset view")
        reset.setToolTip("Look at the target from the selected transducer's side")
        reset.clicked.connect(lambda: self._reset_view(force=True))
        bar.addWidget(reset)
        v.addLayout(bar)

        self._content = QWidget()
        self._content_layout = QVBoxLayout(self._content)
        self._content_layout.setContentsMargins(0, 0, 0, 0)
        self.stack = EmptyStack(self._content)
        v.addWidget(self.stack, 1)
        self.legend = QLabel()
        self.legend.setObjectName("CardSub")
        self.legend.setWordWrap(True)
        v.addWidget(self.legend)

    def _opacity_slider(self, bar, box: QCheckBox, value: float, tip: str) -> QSlider:
        """An opacity slider (10–100 %) right after ``box``, enabled while it is ticked."""
        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setRange(10, 100)
        slider.setValue(round(value * 100))
        slider.setFixedWidth(90)
        slider.setToolTip(tip)
        slider.valueChanged.connect(lambda _v: self._rebuild())
        box.toggled.connect(slider.setEnabled)
        bar.insertWidget(bar.indexOf(box) + 1, slider)
        return slider

    def _load_mask(self) -> None:
        """Read the run's target mask when it is new or changed (one stat per refresh)."""
        path = target_mask_path(self.model)
        try:
            key = (path, path.stat().st_mtime) if path is not None else None
        except OSError:
            key = None
        if key is None:
            self._mask = None
        elif self._mask is None or self._mask[0] != key:
            try:
                vol = pl.Volume(path)
                self._mask = (key, m3.mask_surface(vol.data, vol.affine))
            except Exception:                 # noqa: BLE001 — unreadable mask: not drawn
                self._mask = None
        found = self._mask is not None and len(self._mask[1][1]) > 0
        self.show_mask.setVisible(found)
        self.mask_opacity.setVisible(found)

    # ── refresh ─────────────────────────────────────────────────────────
    def refresh(self) -> None:
        if self.window().isVisible() and not self.isVisible():
            self._stale = True              # another sub-tab is showing: build when shown
            return
        self.v._load_images()
        n = self.model.n_transducers
        self.tx_label.setText(f"Transducer {self.v.ti + 1} of {n} highlighted" if n > 1 else "")
        if self.v.t1 is None:
            self.stack.show_empty("No T1 image yet", "The 3D view shows the head from the subject's "
                                  "segmentation and T1; set them on Setup ▸ Paths.")
            return
        if not can_render():
            self.stack.show_empty("The 3D view needs a display",
                                  "Qt Quick 3D could not be used here (no GPU display, or PyQt6's "
                                  "QtQuick3D module is missing).")
            return
        self.stack.show_content()
        self._load_mask()
        self._ensure_quick()
        self._rebuild()
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
        qml = cache_root() / "placement_3d.qml"          # local disk: QML must load from a file
        qml.parent.mkdir(parents=True, exist_ok=True)
        if not qml.is_file() or qml.read_text(encoding="utf-8") != QML:
            qml.write_text(QML, encoding="utf-8")
        self._quick = QQuickWidget()
        self._quick.setResizeMode(QQuickWidget.ResizeMode.SizeRootObjectToView)
        self._quick.setMinimumHeight(300)
        self._quick.rootContext().setContextProperty("meshes", [])
        self._tip_sink = _TipSink(self._quick)
        self._quick.rootContext().setContextProperty("tipSink", self._tip_sink)
        self._quick.setSource(QUrl.fromLocalFile(str(qml)))
        self._content_layout.addWidget(self._quick, 1)
        from .export import attach_grab, subject_stem
        self.export_button = attach_grab(self._quick, "placement_3d",
                                         subject_stem(lambda: self.model, "placement_3d"))

    # ── scene ───────────────────────────────────────────────────────────
    def meshes(self) -> list[tuple[np.ndarray, np.ndarray, str, float]]:
        """(vertices, faces, colour, opacity) of everything drawn, in RAS mm."""
        out = []
        seg = self.v.seg
        if seg is not None and (self.show_scalp.isChecked() or self.show_brain.isChecked()):
            if self._head is None or self._head[0] is not seg:
                self._head = (seg, m3.head_surfaces(seg))
            surfaces = self._head[1]
            if self.show_brain.isChecked() and len(surfaces["brain"][1]):
                out.append((*surfaces["brain"], BRAIN[0], self.brain_opacity.value() / 100))
            scalp = surfaces["scalp"] if self.show_face.isChecked() else surfaces["scalp_defaced"]
            if self.show_scalp.isChecked() and len(scalp[1]):
                out.append((*scalp, SCALP[0], self.scalp_opacity.value() / 100))
        if self._mask is not None and self.show_mask.isChecked() and len(self._mask[1][1]):
            out.append((*self._mask[1], MASK[0], self.mask_opacity.value() / 100))
        txs = self.v._transducers()
        for i in range(self.model.n_transducers):
            trans, focus = self.v._positions_for(i)
            bowl = self.v._bowl_for(i)
            if trans is None or focus is None or bowl is None or np.linalg.norm(focus - trans) <= 0:
                continue
            selected = i == self.v.ti
            tip = self.model.transducer_tooltip(i)
            axis = (focus - trans) / np.linalg.norm(focus - trans)
            tx = txs[i] if i < len(txs) else {}
            sub = tx.get(tx.get("type")) or {} if tx.get("type") in ("annular", "matrix") else {}
            ids = sub.get("elem_id_mm") if tx.get("type") == "annular" else None
            ods = sub.get("elem_od_mm") if tx.get("type") == "annular" else None
            ids = ids if isinstance(ids, list) else ([ids] if isinstance(ids, (int, float)) else None)
            ods = ods if isinstance(ods, list) else ([ods] if isinstance(ods, (int, float)) else None)
            colours = RINGS_SELECTED if selected else RINGS_OTHER
            for k, (rv, rf) in enumerate(m3.ring_meshes(bowl, ids, ods, trans, axis)):
                out.append(_tipped((rv, rf, colours[k % 2], 1.0 if selected else 0.6), tip))
            if self.show_housing.isChecked():
                housing, known = self.housing_mesh(tx.get("serial"), bowl, trans, axis)
                if selected:
                    self._housing_known = known
                opacity = self.housing_opacity.value() / 100
                if len(housing[1]):
                    out.append(_tipped((*housing, HOUSING if selected else HOUSING_OTHER,
                                        opacity if selected else min(opacity, 0.6)), tip))
            if selected:
                out.append((*m3.exit_plane_disc(bowl, trans, axis), *EXIT_PLANE))
            end = trans + axis * max(float(np.linalg.norm(focus - trans)), bowl.curv_radius
                                     if np.isfinite(bowl.curv_radius) else 0.0)
            out.append(_tipped((*m3.cylinder(trans, end + axis * 10, 0.5), colours[0],
                                1.0 if selected else 0.6), tip))
            out.append((*m3.sphere(focus, 2.5), TARGET if selected else RINGS_OTHER[0], 1.0))
            if selected:                    # a crosshair, so the target stands out inside the brain
                for d in np.eye(3) * 9.0:
                    out.append((*m3.cylinder(focus - d, focus + d, 0.5), TARGET, 1.0))
            if np.isfinite(bowl.curv_radius):
                out.append((*m3.sphere(trans + axis * bowl.curv_radius, 1.5), GEOM_FOCUS, 1.0))
            exp = self.v._expected_bowl_distance(i)
            if exp is not None and selected:
                out.append((*m3.sphere(trans + axis * exp, 2.0), SETTING_FOCUS, 1.0))
        return out

    def housing_mesh(self, serial, bowl, trans, axis):
        """(housing in RAS mm, size on record): a cylinder of the entry's size, else a generic
        cylinder of the aperture, both from the exit plane back past the bowl apex."""
        lib = getattr(self.model, "library", None)
        housing = lib.housing(serial) if lib is not None else {}
        diameter, depth, front, known = m3.housing_size(bowl.aperture, bowl.offset, housing)
        return m3.housing_puck(diameter, depth, front, trans, axis), known

    def _rebuild(self) -> None:
        if self._quick is None:
            return
        # opaque first, transparent last: Qt Quick 3D sorts only within the transparent pass
        meshes = sorted(self.meshes(), key=lambda m: m[3] < 1)
        self._items = [MeshGeometry(v, f, c, a, tip=getattr(m, "tip", "")) for m in meshes
                       for v, f, c, a in [m]]
        self._quick.rootContext().setContextProperty("meshes", self._items)
        root = self._quick.rootObject()
        if root is not None:
            root.setProperty("background", style.active()["surface"])
        mask = self._mask is not None and self.show_mask.isChecked() and len(self._mask[1][1])
        # markers on two lines (positions, transducer parts), the interaction hint on a third;
        # the parenthetical notes go to the tooltip
        self.legend.setText(
            f"<span style='color:{TARGET}'>✛</span> target · "
            + (f"<span style='color:{MASK[0]}'>●</span> target mask · " if mask else "")
            + f"<span style='color:{GEOM_FOCUS}'>●</span> geometric focus · "
            f"<span style='color:{SETTING_FOCUS}'>●</span> focal depth setting<br>"
            f"<span style='color:{RINGS_SELECTED[0]}'>●</span><span style='color:{RINGS_SELECTED[1]}'>●</span> "
            "element rings · "
            + (f"<span style='color:{HOUSING}'>●</span> housing"
               + ("" if self._housing_known else " (generic)")
               + " · " if self.show_housing.isChecked() else "")
            + f"<span style='color:{EXIT_PLANE[0]}'>●</span> exit plane<br>"
            + ("" if self.show_face.isChecked() else "Face scrambled · ")
            + "Drag to turn, scroll to zoom")
        self.legend.setToolTip(
            "Target: the focus position.\nElement rings: the bowl apex sits at the transducer position."
            + ("" if self._housing_known else "\nHousing: generic; the equipment entry gives no size."))

    def _reset_view(self, force: bool = False) -> None:
        """Centre on the selected target, seen from the transducer's side and a little above;
        kept while only positions change, so editing does not move the camera."""
        if self._quick is None or self._quick.rootObject() is None:
            return
        key = (id(self.v.t1), self.v.ti)
        if not force and key == self._view_for:
            return
        self._view_for = key
        trans, focus = self.v.positions()
        if focus is None:
            centre = self.v.t1.to_ras(np.array(self.v.t1.shape) / 2)[0]
            outward = np.array([1.0, 0.6, 0.4])
        else:
            centre = focus
            outward = (trans - focus) if trans is not None and np.linalg.norm(trans - focus) > 0 \
                else np.array([1.0, 0.6, 0.4])
        out = to_scene(outward / np.linalg.norm(outward))[0]
        side = np.cross([0.0, 1.0, 0.0], out)       # oblique, so the bowl's curvature shows
        side = side / np.linalg.norm(side) if np.linalg.norm(side) > 1e-6 else np.zeros(3)
        d = out + 0.6 * side + 0.45 * np.array([0.0, 1.0, 0.0])
        q = look_rotation(d)
        c = to_scene(centre)[0]
        self._quick.rootObject().setProperty(
            "view", [float(c[0]), float(c[1]), float(c[2]), *q, VIEW_DISTANCE_MM])
