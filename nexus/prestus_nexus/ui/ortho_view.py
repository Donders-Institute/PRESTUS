"""Reusable read-only 3-plane orthogonal viewer: an anatomy volume with an optional
colour overlay, both loaded as ``placement.Volume`` and sampled by RAS coordinate
(like PlacementPage._overlay), so the overlay need not share the anatomy's grid or
orientation. Used by the Results tab (acoustic/thermal maps) and the pseudo-CT
viewer — kept separate from PlacementPage, which also draws bowl/transducer geometry
and writes clicks back into the study config.

Volumes are read on a worker thread (the window stays usable while a NIfTI comes off a
network drive) into a small cache, so switching tabs, runs or maps does not read a file
twice; redraws are coalesced and reuse the image artists, so scrolling and the opacity
slider stay responsive."""

from __future__ import annotations

import threading
import time
from collections import OrderedDict
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path

import numpy as np
from PyQt6.QtCore import QObject, Qt, QTimer, pyqtSignal
from PyQt6.QtWidgets import QApplication
from PyQt6.QtWidgets import (QComboBox, QHBoxLayout, QLabel, QSlider, QVBoxLayout, QWidget)

from .. import placement as pl
from .empty_state import EmptyStack
from . import style

try:
    from .mpl_theme import ThemedCanvas as FigureCanvasQTAgg   # drawn in the Nexus theme
    from matplotlib.figure import Figure
    HAVE_MPL = True
except ImportError:          # pragma: no cover
    HAVE_MPL = False

VIEWS = [("Axial", 2, (0, 1)), ("Coronal", 1, (0, 2)), ("Sagittal", 0, (1, 2))]
# Perceptually uniform maps first (magma is the default). No rainbow ('jet': false edges)
# and no diverging map: every overlay shown here (intensity, pressure, MI, temperature, CEM43,
# pseudo-CT) is unsigned, so a map centred on zero would mislead.
COLORMAPS = ["magma", "inferno", "viridis", "plasma", "hot"]

def _fmt(v: float | None) -> str:
    """A map value for the readout: 3 significant digits, no exponent below 1e5."""
    if v is None:
        return ""
    return f"{v:.3g}" if abs(v) < 1e5 else f"{v:.3e}"


# ── volume cache ─────────────────────────────────────────────────────────
_CACHE_SIZE = 6              # a T1 plus a few maps; each holds its data in memory
_cache: "OrderedDict[tuple, pl.Volume]" = OrderedDict()
_inflight: dict[tuple, Future] = {}
_lock = threading.Lock()
_executor = ThreadPoolExecutor(max_workers=2, thread_name_prefix="nexus-volume")


def _key(path: str | Path) -> tuple:
    p = Path(path).resolve()
    st = p.stat()
    return (str(p), st.st_mtime_ns, st.st_size)


def _load(key: tuple) -> pl.Volume:
    try:
        vol = pl.Volume(key[0])
        vol.canon                         # reads the voxels (and reorients them)
        with _lock:
            for k in [k for k in _cache if k[0] == key[0]]:   # an older version of the file
                del _cache[k]
            _cache[key] = vol
            while len(_cache) > _CACHE_SIZE:
                _cache.popitem(last=False)
        return vol
    finally:
        with _lock:
            _inflight.pop(key, None)


def load_volume_async(path: str | Path) -> Future:
    """A Future of the ``placement.Volume`` for ``path``, with its voxels read. Done at
    once when the cache holds the file unchanged (same size and modification time);
    else read on a worker thread, once however many callers ask."""
    key = _key(path)
    with _lock:
        vol = _cache.get(key)
        if vol is not None:
            _cache.move_to_end(key)
            fut: Future = Future()
            fut.set_result(vol)
            return fut
        fut = _inflight.get(key)
        if fut is None:
            fut = _inflight[key] = _executor.submit(_load, key)
        return fut


def load_volume(path: str | Path) -> pl.Volume:
    """``load_volume_async(path)``, waited for."""
    return load_volume_async(path).result()


def cached_volume(path: str | Path) -> pl.Volume | None:
    """The cached volume for ``path`` if it is loaded and unchanged, else None."""
    try:
        key = _key(path)
    except OSError:
        return None
    with _lock:
        return _cache.get(key)


class _Relay(QObject):
    """Hands a finished load from the worker thread to the viewer's (GUI) thread."""
    done = pyqtSignal(object, object)     # (role, token), Future


class OrthoViewer(QWidget):
    """Axial/coronal/sagittal slices through an anatomy volume, with an optional
    colour overlay (own colormap + opacity). No click-to-place; view only."""

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.anat: pl.Volume | None = None
        self.overlay: pl.Volume | None = None
        self.cursor = np.zeros(3)
        self.markers: list[tuple[np.ndarray, str, str]] = []   # (RAS mm, marker, colour)
        self.marker_texts: list[str] = []         # hover text per marker ('' = none)
        self._anat_lim: tuple[float, float] | None = None
        self._over_lim: tuple[float, float] | None = None
        self._over_same_grid = False
        self._over_slices: dict[tuple[int, int], np.ndarray] = {}
        self._artists: list[dict] | None = None   # per view: anat/overlay images, markers
        self._peak: dict[int, np.ndarray | None] = {}   # id(overlay) -> RAS of its maximum
        self._pending: dict[str, object] = {}     # role ('anat'/'overlay') -> load token
        self._when_loaded: list = []              # callbacks once nothing is pending
        self.unit = ""                            # overlay unit, for the colour bar and readout
        self.zero_based = False                   # colour scale 0..max instead of 1..99th percentile
        self._colorbar = None
        self._relay = _Relay(self)
        self._relay.done.connect(self._loaded, Qt.ConnectionType.QueuedConnection)

        # What the empty state says when there is no anatomy: (title, body, action, callback);
        # the host page sets it (set_empty_message), so the panel names the missing input.
        self._empty_msg: tuple = ("No image to show", "", None, None)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        content = QWidget()
        self.stack = EmptyStack(content)
        outer.addWidget(self.stack)
        v = QVBoxLayout(content)
        v.setContentsMargins(4, 4, 4, 4)
        v.setSpacing(4)

        tb = QHBoxLayout()
        # overlay controls: shown only while there is an overlay to adjust
        self.overlay_controls = QWidget()
        oc = QHBoxLayout(self.overlay_controls)
        oc.setContentsMargins(0, 0, 0, 0)
        tb.addWidget(self.overlay_controls)
        oc.addWidget(QLabel("Overlay colormap"))
        self.cmap = QComboBox()
        self.cmap.addItems(COLORMAPS)
        self.cmap.setFixedWidth(110)
        self.cmap.currentTextChanged.connect(self._request_draw)
        oc.addWidget(self.cmap)
        oc.addWidget(QLabel("Opacity"))
        self.alpha = QSlider(Qt.Orientation.Horizontal, minimum=0, maximum=100, value=50)
        self.alpha.setMaximumWidth(140)
        self.alpha.valueChanged.connect(self._request_draw)
        oc.addWidget(self.alpha)
        oc.addWidget(QLabel("Hide below"))
        self.threshold = QSlider(Qt.Orientation.Horizontal, minimum=0, maximum=90, value=0,
                                 singleStep=5, pageStep=10)
        self.threshold.setMaximumWidth(120)
        self.threshold.setToolTip("Hide overlay values below this share of the colour scale's top")
        self.threshold.valueChanged.connect(self._request_draw)
        oc.addWidget(self.threshold)
        self.threshold_label = QLabel()
        self.threshold_label.setMinimumWidth(110)
        oc.addWidget(self.threshold_label)
        self.overlay_controls.hide()
        tb.addStretch(1)
        self.status = QLabel()
        self.status.setObjectName("CardSub")
        tb.addWidget(self.status)
        v.addLayout(tb)

        # Several setters in a row (anatomy, markers, overlay, centre) draw once.
        self._draw_timer = QTimer(self, singleShot=True, interval=0)
        self._draw_timer.timeout.connect(self._draw)

        if HAVE_MPL:
            self.fig = Figure(figsize=(9, 3.2))
            self.canvas = FigureCanvasQTAgg(self.fig, export_name="slices")
            self.axes = [self.fig.add_subplot(1, 3, i + 1) for i in range(3)]
            # fixed margins (not tight_layout): room for the colour bar below the slices
            self.fig.subplots_adjust(left=0.02, right=0.98, top=0.93, bottom=0.14, wspace=0.06)
            self.canvas.mpl_connect("scroll_event", self._scrolled)
            self.canvas.mpl_connect("button_press_event", self._clicked)
            from .tooltip import MarkerTips
            self.marker_tips = MarkerTips(self.canvas)     # hover a marker: its transducer's name
            v.addWidget(self.canvas, 1)
            self.readout = QLabel()
            self.readout.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            v.addWidget(self.readout)
            self.stack.show_empty(*self._empty_msg)       # never blank axes with 0–1 ticks
        else:
            self._empty_msg = ("This viewer needs matplotlib and nibabel",
                               "Install them in the Nexus's Python environment and restart.", None, None)
            self.stack.show_empty(*self._empty_msg)

    def set_empty_message(self, title: str, body: str = "", action: str | None = None,
                          callback=None) -> None:
        """What to say (and offer) while there is no anatomy to show."""
        if HAVE_MPL:
            self._empty_msg = (title, body, action, callback)
            if self.anat is None and not self._pending:
                self.stack.show_empty(*self._empty_msg)

    # ── loading ─────────────────────────────────────────────────────────
    def set_anatomy(self, path: str | Path | None) -> None:
        """Show this anatomy (loaded in the background unless cached)."""
        self._request("anat", path)

    def set_overlay(self, path: str | Path | None, unit: str | None = None,
                    zero_based: bool | None = None, threshold_pct: int | None = None) -> None:
        """Show this colour overlay (loaded in the background unless cached). ``unit`` labels
        the colour bar and readout; ``zero_based`` scales colours from 0 to the map's
        maximum (intensity, pressure: the peak is not clipped) instead of its 1st..99th
        percentile; ``threshold_pct`` sets 'Hide below'."""
        if unit is not None:
            self.unit = unit
        if zero_based is not None and zero_based != self.zero_based:
            self.zero_based = zero_based
            self.overlay = None                   # recompute the colour limits
        if threshold_pct is not None:
            self.threshold.setValue(int(threshold_pct))
        self._request("overlay", path)

    def when_loaded(self, callback) -> None:
        """Call ``callback()`` once the requested anatomy and overlay are shown (now,
        if nothing is loading)."""
        if self._pending:
            self._when_loaded.append(callback)
        else:
            callback()

    def is_loading(self) -> bool:
        return bool(self._pending)

    def wait(self, timeout: float = 120.0) -> None:
        """Block (processing events) until pending loads are shown; for tests/scripts."""
        end = time.monotonic() + timeout
        while self._pending and time.monotonic() < end:
            QApplication.processEvents()
            time.sleep(0.005)

    def _request(self, role: str, path: str | Path | None) -> None:
        token = object()
        self._pending[role] = token
        fut = None
        if path is not None and Path(path).is_file():
            try:
                fut = load_volume_async(path)
            except OSError as exc:
                self.status.setText(f"Cannot read {Path(path).name}: {exc}")
        if fut is None:
            self._finish(role, token, None)
        elif fut.done():
            self._finish(role, token, fut)
        else:
            self.status.setText(f"Loading {Path(path).name} …")
            relay = self._relay

            def done(f, r=role, t=token):
                try:
                    relay.done.emit((r, t), f)
                except RuntimeError:  # best-effort: viewer closed meanwhile
                    pass
            fut.add_done_callback(done)

    def _loaded(self, role_token, fut: Future) -> None:
        role, token = role_token
        if self._pending.get(role) is token:      # else superseded by a newer request
            self._finish(role, token, fut)

    def _finish(self, role: str, token, fut: Future | None) -> None:
        vol = None
        if fut is not None:
            try:
                vol = fut.result()
                self.status.setText("")
            except Exception as exc:  # reported: in the status line
                self.status.setText(f"Cannot read {'anatomy' if role == 'anat' else 'overlay'}: {exc}")
        if self._pending.get(role) is token:
            del self._pending[role]
        if self._pending and vol is not None:
            self.status.setText("Loading …")
        (self._apply_anatomy if role == "anat" else self._apply_overlay)(vol)
        if not self._pending:
            callbacks, self._when_loaded = self._when_loaded, []
            for cb in callbacks:
                cb()

    def _apply_anatomy(self, vol: pl.Volume | None) -> None:
        if vol is not self.anat:
            self.anat = vol
            self._artists = None
            self._over_slices.clear()
            self._anat_lim = None
            if vol is not None:
                self.cursor = vol.to_ras(np.array(vol.shape) / 2)[0]
                lo, hi = np.percentile(vol.canon[::4, ::4, ::4], [1, 99.5])
                self._anat_lim = (float(lo), float(hi))
            self._match_grids()
        self._request_draw()

    def _apply_overlay(self, vol: pl.Volume | None) -> None:
        if vol is not self.overlay:
            self.overlay = vol
            self._over_slices.clear()
            self._over_lim = None
            if vol is not None:
                d = vol.data
                nz = d[(d != 0) & np.isfinite(d)]
                if nz.size and self.zero_based:
                    self._over_lim = (0.0, float(nz.max()))
                elif nz.size:
                    lo, hi = np.percentile(nz, [1, 99])
                    self._over_lim = (float(lo), float(hi))
            self._match_grids()
        self._request_draw()

    def _match_grids(self) -> None:
        """Overlay on the anatomy's grid (the usual T1w maps): slice it directly
        instead of resampling it by RAS coordinate."""
        a, o = self.anat, self.overlay
        self._over_same_grid = (a is not None and o is not None
                                and tuple(a.canon.shape[:3]) == tuple(o.canon.shape[:3])
                                and np.allclose(a.canon_affine, o.canon_affine, atol=1e-4))

    def set_markers(self, markers: list[tuple]) -> None:
        """Points to draw on every view: (RAS mm, matplotlib marker, colour[, hover text]).
        A point within two voxels of a slice is drawn solid, otherwise faint (projected)."""
        self.markers = [(np.asarray(mk[0], float), mk[1], mk[2]) for mk in markers]
        self.marker_texts = [str(mk[3]) if len(mk) > 3 else "" for mk in markers]
        self._request_draw()

    def centre_on(self, ras) -> None:
        """Move the crosshair (and so the three slices) to a RAS point."""
        if ras is not None:
            self.cursor = np.asarray(ras, float)
            self._request_draw()

    def overlay_peak(self) -> np.ndarray | None:
        """RAS of the overlay's maximum, or None."""
        if self.overlay is None:
            return None
        key = id(self.overlay)
        if key not in self._peak:
            self._peak = {key: self._find_peak()}
        return self._peak[key]

    def _find_peak(self) -> np.ndarray | None:
        d = np.nan_to_num(np.asarray(self.overlay.data, dtype=float))
        if not d.size or not np.any(d):
            return None
        return self.overlay.to_ras(np.unravel_index(int(np.argmax(d)), d.shape))[0]

    # ── drawing ─────────────────────────────────────────────────────────
    def _request_draw(self, *_args) -> None:
        self._draw_timer.start()

    def _overlay_slice(self, plane: int, k: int, shape2d: tuple[int, int]) -> np.ndarray | None:
        if self.overlay is None or self.anat is None:
            return None
        cached = self._over_slices.get((plane, k))
        if cached is not None:
            return cached
        if self._over_same_grid:
            values = np.asarray(np.take(self.overlay.canon, k, axis=plane), dtype=float)
        else:
            ii, jj = np.meshgrid(np.arange(shape2d[0]), np.arange(shape2d[1]), indexing="ij")
            vox = np.zeros((ii.size, 3))
            others = [a for a in range(3) if a != plane]
            vox[:, others[0]] = ii.ravel()
            vox[:, others[1]] = jj.ravel()
            vox[:, plane] = k
            values = self.overlay.sample(self.anat.canon_to_ras(vox)).astype(float).reshape(shape2d)
        masked = np.where(values != 0, values, np.nan)
        if len(self._over_slices) > 48:
            self._over_slices.clear()
        self._over_slices[(plane, k)] = masked
        return masked

    def _draw(self) -> None:
        self._draw_timer.stop()
        if not HAVE_MPL:
            return
        self.overlay_controls.setVisible(self.overlay is not None)
        if self.anat is None:
            self._artists = None
            for ax in self.axes:
                ax.clear()
                ax.set_xticks([])
                ax.set_yticks([])
            if self._pending:
                self.stack.show_empty("Loading the image …", "Large NIfTI files on a network "
                                      "drive can take a few seconds.")
            else:
                self.stack.show_empty(*self._empty_msg)
            return
        self.stack.show_content()
        data = self.anat.canon
        vs = np.sqrt((self.anat.canon_affine[:3, :3] ** 2).sum(axis=0))
        c = self.anat.ras_to_canon(self.cursor)[0]
        alpha = self.alpha.value() / 100
        show_over = alpha > 0 and self._over_lim is not None and self._over_lim[1] > self._over_lim[0]
        thr = None
        if show_over and self.threshold.value() > 0:
            lo, hi = self._over_lim
            thr = lo + (hi - lo) * self.threshold.value() / 100
        self.threshold_label.setText(f"≥ {_fmt(thr)} {self.unit}".rstrip() if thr is not None else "")
        if self._artists is None:
            self._artists = []
            for ax, (title, plane, (h, w)) in zip(self.axes, VIEWS):
                ax.clear()
                ax.set_xticks([])
                ax.set_yticks([])
                shape2d = tuple(s for i, s in enumerate(data.shape[:3]) if i != plane)
                blank = np.zeros(shape2d).T
                anat_im = ax.imshow(blank, cmap="gray", origin="lower", vmin=self._anat_lim[0],
                                    vmax=self._anat_lim[1], aspect=vs[w] / vs[h],
                                    interpolation="nearest")
                over_im = ax.imshow(np.full_like(blank, np.nan, dtype=float), origin="lower",
                                    aspect=vs[w] / vs[h], interpolation="nearest")
                hline = ax.axhline(0, color="#4fa3e0", lw=0.6, alpha=0.6)
                vline = ax.axvline(0, color="#4fa3e0", lw=0.6, alpha=0.6)
                ax.set_title(title, fontsize=style.PLOT_PT["title"])
                ax.set_xlim(-0.5, shape2d[0] - 0.5)
                ax.set_ylim(-0.5, shape2d[1] - 0.5)
                self._artists.append({"anat": anat_im, "over": over_im, "k": None, "marks": [],
                                      "cross": (hline, vline)})
        for ax, art, (title, plane, (h, w)) in zip(self.axes, self._artists, VIEWS):
            k = int(np.clip(round(c[plane]), 0, data.shape[plane] - 1))
            if art["k"] != k:
                sl = np.take(data, k, axis=plane)
                art["anat"].set_data(sl.T)
                art["k"] = k
                art["shape"] = sl.shape
            over = self._overlay_slice(plane, k, art["shape"]) if show_over else None
            if over is not None:
                if thr is not None:
                    with np.errstate(invalid="ignore"):
                        over = np.where(over >= thr, over, np.nan)
                art["over"].set_data(over.T)
                art["over"].set_cmap(self.cmap.currentText())
                art["over"].set_clim(*self._over_lim)
                art["over"].set_alpha(alpha)
                art["over"].set_visible(True)
            else:
                art["over"].set_visible(False)
            art["cross"][0].set_ydata([c[w], c[w]])
            art["cross"][1].set_xdata([c[h], c[h]])
            for line in art["marks"]:
                line.remove()
            art["marks"] = []
            self.marker_tips.clear(ax)
            for (ras, marker, colour), tip in zip(self.markers, self.marker_texts):
                m = self.anat.ras_to_canon(ras)[0]
                near = abs(m[plane] - k) <= 2
                art["marks"] += ax.plot(m[h], m[w], marker, color=colour, markersize=8 if near else 6,
                                        markeredgewidth=1.5, alpha=1.0 if near else 0.35,
                                        linestyle="none")
                self.marker_tips.add(ax, m[h], m[w], tip)
        self._update_colorbar(show_over)
        self._update_readout()
        self.canvas.draw_idle()

    def _update_colorbar(self, show: bool) -> None:
        """One horizontal colour bar under the slices, with the overlay's unit."""
        if not show:
            if self._colorbar is not None:
                self._colorbar.ax.set_visible(False)
            return
        mappable = self._artists[0]["over"]
        if self._colorbar is None:
            cax = self.fig.add_axes([0.35, 0.04, 0.3, 0.03])
            self._colorbar = self.fig.colorbar(mappable, cax=cax, orientation="horizontal")
            self._colorbar.ax.tick_params(labelsize=style.PLOT_PT["small"])
        else:
            self._colorbar.update_normal(mappable)
        self._colorbar.ax.set_visible(True)
        self._colorbar.set_alpha(1.0)
        self._colorbar.set_label(f"[{self.unit}]" if self.unit else "", fontsize=style.PLOT_PT["normal"])

    def _update_readout(self) -> None:
        """RAS of the crosshair and the overlay's value there."""
        if self.anat is None:
            self.readout.setText("")
            return
        r = self.cursor
        text = f"Crosshair RAS [{r[0]:.1f}, {r[1]:.1f}, {r[2]:.1f}] mm"
        if self.overlay is not None:
            value = float(np.asarray(self.overlay.sample(r[None]), dtype=float).ravel()[0])
            text += (f" · overlay {_fmt(value)} {self.unit}".rstrip() if np.isfinite(value)
                     else " · overlay: outside the map")
        text += "  (click a slice to move the crosshair)"
        self.readout.setText(text)

    def _clicked(self, event) -> None:
        """A click moves the crosshair to that point (view only; nothing is written)."""
        if self.anat is None or event.inaxes not in self.axes or event.xdata is None:
            return
        _, plane, (h, w) = VIEWS[self.axes.index(event.inaxes)]
        c = self.anat.ras_to_canon(self.cursor)[0]
        c[h], c[w] = event.xdata, event.ydata
        self.cursor = self.anat.canon_to_ras(c[None])[0]
        self._request_draw()

    def _scrolled(self, event) -> None:
        if self.anat is None or event.inaxes not in self.axes:
            return
        _, plane, _ = VIEWS[self.axes.index(event.inaxes)]
        c = self.anat.ras_to_canon(self.cursor)[0]
        c[plane] += 1 if event.button == "up" else -1
        self.cursor = self.anat.canon_to_ras(c[None])[0]
        self._request_draw()
