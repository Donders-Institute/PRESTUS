"""Figure export: a small download button pinned to the top-right corner of a figure.

Three kinds of figure, one button (``ExportButton``):

- matplotlib (``MplExporter``): SVG and PDF (vector), PNG at 300 or 600 dpi, copy. ThemedCanvas
  (mpl_theme.py) attaches the button itself.
- hand-painted QWidgets (``PaintedExportable``): the widget paints through
  ``paint(qp, rect, palette)``, so the same code draws to the screen, an SVG, a PDF and a PNG
  at 3x or 6x.
- 3D views (``GrabExporter``): PNG of the rendered frame at the screen's pixel ratio, copy.

Unless "Keep current theme" is ticked (remembered for the session), figures are exported in
the light palette on white, for papers and slides, and the screen keeps its theme. The last
folder used is remembered for the session.
"""

from __future__ import annotations

import io
from pathlib import Path
from typing import Callable

from PyQt6.QtCore import QEvent, QMarginsF, QObject, QPoint, QRect, QRectF, QSize, QSizeF, Qt
from PyQt6.QtGui import (QAction, QColor, QGuiApplication, QImage, QPageLayout, QPageSize, QPainter,
                         QPdfWriter)
from PyQt6.QtWidgets import QFileDialog, QMenu, QMessageBox, QToolButton, QWidget

from . import style

# Session state shared by every button.
_SESSION = {"dir": None, "keep_theme": False}

# (menu text, format, dpi or scale)
VECTOR = [("SVG (vector)…", "svg", None), ("PDF (vector)…", "pdf", None)]
RASTER = [("PNG, 300 dpi…", "png", 300), ("PNG, 600 dpi…", "png", 600)]


def keep_theme() -> bool:
    return bool(_SESSION["keep_theme"])


def set_keep_theme(on: bool) -> None:
    _SESSION["keep_theme"] = bool(on)


def export_palette() -> dict:
    """The palette figures are exported in: light, unless the user keeps the current theme."""
    return style.active() if keep_theme() else style.C


def subject_stem(model, view: str) -> Callable[[], str]:
    """A filename hint ``sub-XXX_<view>[_desc-<affix>]`` from the model's subject and output
    affix, or just ``<view>`` when there is no subject. ``model`` may be a zero-argument
    callable returning the model, for pages whose model is replaced when a config opens."""
    def hint() -> str:
        from ..study import sanitize_affix
        try:
            m = model() if callable(model) and not hasattr(model, "get") else model
            sid = m.get("subject_id")
            affix = sanitize_affix(m.get("io.output_affix") or "")
        except Exception:  # best-effort: a file name hint must never fail
            return view
        if isinstance(sid, (int, float)) and not isinstance(sid, bool):
            return f"sub-{int(sid):03d}_{view}{affix}"
        return f"{view}{affix}"
    return hint


# ── exporters ───────────────────────────────────────────────────────────
class MplExporter:
    """Writes a matplotlib figure; re-themes it to the export palette for the file only."""

    formats = VECTOR + RASTER

    def __init__(self, canvas):
        self.canvas = canvas

    def _themed(self, write: Callable[[dict], None]) -> None:
        from . import mpl_theme
        fig = self.canvas.figure
        pal = export_palette()
        try:
            mpl_theme.apply(fig, pal)
            write(pal)
        finally:
            try:
                mpl_theme.apply(fig, style.active())
            finally:
                self.canvas.draw_idle()

    def save(self, path: Path, fmt: str, dpi: int | None) -> None:
        def write(pal):
            self.canvas.figure.savefig(str(path), format=fmt, dpi=dpi or 300, bbox_inches="tight",
                                       facecolor=pal["surface"])
        self._themed(write)

    def image(self) -> QImage:
        buf = io.BytesIO()

        def write(pal):
            self.canvas.figure.savefig(buf, format="png", dpi=200, bbox_inches="tight",
                                       facecolor=pal["surface"])
        self._themed(write)
        return QImage.fromData(buf.getvalue(), "PNG")


class PaintedExportable:
    """Mixin for a QWidget drawn by ``paint(qp, rect, palette)``: paintEvent paints the screen
    in the active palette, and the export paints the same picture to SVG, PDF or PNG."""

    export_name = "figure"

    def paint(self, qp: QPainter, rect: QRectF, c: dict) -> None:      # pragma: no cover
        raise NotImplementedError

    def paintEvent(self, _event) -> None:  # noqa: N802
        qp = QPainter(self)
        qp.setRenderHint(QPainter.RenderHint.Antialiasing)
        try:
            self.paint(qp, QRectF(self.rect()), style.active())
        finally:
            qp.end()


def _have_qtsvg() -> bool:
    try:
        from PyQt6.QtSvg import QSvgGenerator  # noqa: F401
        return True
    except ImportError:                     # pragma: no cover — PyQt6 without QtSvg: PDF only
        return False


class PaintedExporter:
    formats = ([f for f in VECTOR if f[1] != "svg" or _have_qtsvg()]
               + [("PNG, 3× screen…", "png", 3), ("PNG, 6× screen…", "png", 6)])

    def __init__(self, widget):
        self.widget = widget

    def _size(self) -> QSize:
        s = self.widget.size()
        return s if s.width() > 0 and s.height() > 0 else self.widget.sizeHint()

    def _paint(self, device, scale: float = 1.0) -> None:
        size = self._size()
        qp = QPainter(device)
        try:
            qp.setRenderHint(QPainter.RenderHint.Antialiasing)
            if scale != 1.0:
                qp.scale(scale, scale)
            self.widget.paint(qp, QRectF(0, 0, size.width(), size.height()), export_palette())
        finally:
            qp.end()

    def raster(self, scale: float) -> QImage:
        size = self._size()
        img = QImage(round(size.width() * scale), round(size.height() * scale),
                     QImage.Format.Format_ARGB32)
        img.fill(QColor(export_palette()["surface"]))
        self._paint(img, scale)
        return img

    def save(self, path: Path, fmt: str, dpi: int | None) -> None:
        size = self._size()
        if fmt == "svg":
            from PyQt6.QtSvg import QSvgGenerator
            gen = QSvgGenerator()
            gen.setFileName(str(path))
            gen.setSize(size)
            gen.setViewBox(QRect(0, 0, size.width(), size.height()))
            gen.setTitle(getattr(self.widget, "export_name", "figure"))
            self._paint(gen)
        elif fmt == "pdf":
            pdf = QPdfWriter(str(path))
            pdf.setResolution(72)                       # 1 widget px = 1 pt
            pdf.setPageSize(QPageSize(QSizeF(size.width(), size.height()), QPageSize.Unit.Point))
            pdf.setPageMargins(QMarginsF(0, 0, 0, 0), QPageLayout.Unit.Point)
            self._paint(pdf)
        else:
            img = self.raster(float(dpi or 3))
            img.setDotsPerMeterX(round(96 * (dpi or 3) / 0.0254))
            img.setDotsPerMeterY(round(96 * (dpi or 3) / 0.0254))
            if not img.save(str(path), "PNG"):
                raise OSError(f"could not write {path}")

    def image(self) -> QImage:
        return self.raster(3.0)


class GrabExporter:
    """3D views: the rendered frame as PNG, at the screen's device pixel ratio."""

    formats = [("PNG (as shown)…", "png", None)]

    def __init__(self, widget, grab: Callable[[], QImage] | None = None):
        self.widget = widget
        self.grab = grab

    def image(self) -> QImage:
        if self.grab is not None:
            img = self.grab()
        elif hasattr(self.widget, "grabFramebuffer"):
            img = self.widget.grabFramebuffer()
        else:
            img = self.widget.grab().toImage()
        return img

    def save(self, path: Path, fmt: str, dpi: int | None) -> None:
        if not self.image().save(str(path), "PNG"):
            raise OSError(f"could not write {path}")


# ── the button ──────────────────────────────────────────────────────────
class ExportButton(QToolButton):
    """A small download icon pinned to the top-right corner of ``owner``, with the export menu.

    ``exporter`` is an MplExporter, PaintedExporter or GrabExporter; ``view`` names the figure
    (the default filename); ``filename_hint`` (callable -> filename stem) overrides it, e.g.
    ``subject_stem(model, "beam")``."""

    MARGIN = 4

    def __init__(self, owner: QWidget, exporter, view: str = "figure",
                 filename_hint: Callable[[], str] | None = None):
        super().__init__(owner)
        self.owner = owner
        self.exporter = exporter
        self.view = view
        self.filename_hint = filename_hint
        self._available = True
        self.setObjectName("ExportButton")
        self.setToolTip("Export figure")
        self.setAutoRaise(True)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.setFixedSize(24, 24)
        self.setIconSize(QSize(16, 16))
        self._menu = QMenu(self)
        self._menu.aboutToShow.connect(self._fill_menu)
        self.setMenu(self._menu)
        self.restyle()
        owner.installEventFilter(self)
        self._place()
        self.raise_()

    # appearance and position
    def restyle(self) -> None:
        from . import icons
        c = style.active()
        self.setIcon(icons.icon("download", c["ink3"]))
        # a plain tinted box behind the glyph so it reads over plots and images
        self.setStyleSheet(
            f"QToolButton#ExportButton {{ background: {_rgba(c['surface'], 0.85)}; border: 1px solid {c['line']};"
            f" border-radius: 5px; padding: 0; }}"
            f"QToolButton#ExportButton:hover {{ background: {c['surface3']}; }}"
            "QToolButton#ExportButton::menu-indicator { image: none; width: 0; }")

    def _place(self) -> None:
        w = self.owner.width()
        self.move(QPoint(max(0, w - self.width() - self.MARGIN), self.MARGIN))
        self.raise_()

    def eventFilter(self, obj: QObject, ev: QEvent) -> bool:  # noqa: N802
        if obj is self.owner:
            t = ev.type()
            if t in (QEvent.Type.Resize, QEvent.Type.Show):
                self._place()
            elif t in (QEvent.Type.PaletteChange, QEvent.Type.StyleChange):
                self.restyle()
        return False

    def set_available(self, on: bool) -> None:
        """Hide the button while the figure has nothing to export (empty state, no data)."""
        on = bool(on)
        if on == self._available and self.isHidden() != on:
            return                              # unchanged (called from paint and draw)
        self._available = on
        self.setVisible(on)
        if on:
            self._place()

    # menu
    def _fill_menu(self) -> None:
        self.restyle()                                  # follows a theme switch
        m = self._menu
        m.clear()
        for text, fmt, dpi in self.exporter.formats:
            act = m.addAction(text)
            act.triggered.connect(lambda _c=False, f=fmt, d=dpi: self.export_to_dialog(f, d))
        m.addSeparator()
        m.addAction("Copy image").triggered.connect(lambda _c=False: self.copy())
        if not isinstance(self.exporter, GrabExporter):
            m.addSeparator()
            keep = QAction("Keep current theme", m)
            keep.setCheckable(True)
            keep.setChecked(keep_theme())
            keep.setToolTip("Export in the theme on screen; unticked, figures are exported "
                            "in the light theme on white")
            keep.toggled.connect(set_keep_theme)
            m.addAction(keep)

    def default_name(self, fmt: str, dpi: int | None = None) -> str:
        stem = self.view
        if self.filename_hint is not None:
            try:
                stem = self.filename_hint() or self.view
            except Exception:  # best-effort: a file name hint must never fail
                stem = self.view
        return f"{stem}.{fmt}"

    def export_to_dialog(self, fmt: str, dpi: int | None) -> Path | None:
        folder = _SESSION["dir"] or str(Path.home())
        start = str(Path(folder) / self.default_name(fmt, dpi))
        flt = {"svg": "SVG (*.svg)", "pdf": "PDF (*.pdf)", "png": "PNG (*.png)"}[fmt]
        path, _ = QFileDialog.getSaveFileName(self.window(), "Export figure", start, flt)
        if not path:
            return None
        path = Path(path)
        if path.suffix.lower() != f".{fmt}":
            path = path.with_suffix(f".{fmt}")
        try:
            self.export(path, fmt, dpi)
        except Exception as exc:  # reported: in a dialog; exporting must never crash the app
            QMessageBox.warning(self.window(), "Export figure", f"Could not export the figure:\n{exc}")
            return None
        return path

    def export(self, path: Path, fmt: str | None = None, dpi: int | None = None) -> Path:
        """Write the figure to ``path`` (format from the suffix unless given)."""
        path = Path(path)
        fmt = (fmt or path.suffix.lstrip(".")).lower()
        self.exporter.save(path, fmt, dpi)
        _SESSION["dir"] = str(path.parent)
        return path

    def copy(self) -> None:
        img = self.exporter.image()
        cb = QGuiApplication.clipboard()
        if cb is not None and not img.isNull():
            cb.setImage(img)


def _rgba(hex_colour: str, alpha: float) -> str:
    q = QColor(hex_colour)
    return f"rgba({q.red()}, {q.green()}, {q.blue()}, {round(alpha * 255)})"


def attach_painted(widget: QWidget, view: str | None = None,
                   filename_hint: Callable[[], str] | None = None) -> ExportButton:
    """The export button on a PaintedExportable widget."""
    return ExportButton(widget, PaintedExporter(widget), view or getattr(widget, "export_name", "figure"),
                        filename_hint)


def attach_grab(widget: QWidget, view: str, filename_hint: Callable[[], str] | None = None,
                grab: Callable[[], QImage] | None = None) -> ExportButton:
    """The export button on a 3D (or any rendered) widget: PNG and copy only."""
    return ExportButton(widget, GrabExporter(widget, grab), view, filename_hint)
