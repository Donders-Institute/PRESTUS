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

try:
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
    from matplotlib.figure import Figure
    HAVE_MPL = True
except ImportError:          # pragma: no cover
    HAVE_MPL = False

VIEWS = [("Axial", 2, (0, 1)), ("Coronal", 1, (0, 2)), ("Sagittal", 0, (1, 2))]
COLORMAPS = ["hot", "plasma", "viridis", "jet", "coolwarm"]

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
        self._anat_lim: tuple[float, float] | None = None
        self._over_lim: tuple[float, float] | None = None
        self._over_same_grid = False
        self._over_slices: dict[tuple[int, int], np.ndarray] = {}
        self._artists: list[dict] | None = None   # per view: anat/overlay images, markers
        self._peak: dict[int, np.ndarray | None] = {}   # id(overlay) -> RAS of its maximum
        self._pending: dict[str, object] = {}     # role ('anat'/'overlay') -> load token
        self._when_loaded: list = []              # callbacks once nothing is pending
        self._relay = _Relay(self)
        self._relay.done.connect(self._loaded, Qt.ConnectionType.QueuedConnection)

        v = QVBoxLayout(self)
        v.setContentsMargins(4, 4, 4, 4)
        v.setSpacing(4)

        tb = QHBoxLayout()
        tb.addWidget(QLabel("Overlay colormap"))
        self.cmap = QComboBox()
        self.cmap.addItems(COLORMAPS)
        self.cmap.setFixedWidth(110)
        self.cmap.currentTextChanged.connect(self._request_draw)
        tb.addWidget(self.cmap)
        tb.addWidget(QLabel("Opacity"))
        self.alpha = QSlider(Qt.Orientation.Horizontal, minimum=0, maximum=100, value=50)
        self.alpha.setMaximumWidth(140)
        self.alpha.valueChanged.connect(self._request_draw)
        tb.addWidget(self.alpha)
        tb.addStretch(1)
        self.status = QLabel()
        self.status.setStyleSheet("color: gray;")
        tb.addWidget(self.status)
        v.addLayout(tb)

        # Several setters in a row (anatomy, markers, overlay, centre) draw once.
        self._draw_timer = QTimer(self, singleShot=True, interval=0)
        self._draw_timer.timeout.connect(self._draw)

        if HAVE_MPL:
            self.fig = Figure(figsize=(9, 3.2), tight_layout=True)
            self.canvas = FigureCanvasQTAgg(self.fig)
            self.axes = [self.fig.add_subplot(1, 3, i + 1) for i in range(3)]
            self.canvas.mpl_connect("scroll_event", self._scrolled)
            v.addWidget(self.canvas, 1)
        else:
            v.addWidget(QLabel("Install matplotlib and nibabel to use this viewer."), 1)

    # ── loading ─────────────────────────────────────────────────────────
    def set_anatomy(self, path: str | Path | None) -> None:
        """Show this anatomy (loaded in the background unless cached)."""
        self._request("anat", path)

    def set_overlay(self, path: str | Path | None) -> None:
        """Show this colour overlay (loaded in the background unless cached)."""
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
                except RuntimeError:          # viewer closed meanwhile
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
            except Exception as exc:
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
                if nz.size:
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

    def set_markers(self, markers: list[tuple[np.ndarray, str, str]]) -> None:
        """Points to draw on every view: (RAS mm, matplotlib marker, colour). A point
        within two voxels of a slice is drawn solid, otherwise faint (projected)."""
        self.markers = [(np.asarray(r, float), m, c) for r, m, c in markers]
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
        if self.anat is None:
            self._artists = None
            for ax in self.axes:
                ax.clear()
                ax.set_xticks([])
                ax.set_yticks([])
            self.canvas.draw_idle()
            return
        data = self.anat.canon
        vs = np.sqrt((self.anat.canon_affine[:3, :3] ** 2).sum(axis=0))
        c = self.anat.ras_to_canon(self.cursor)[0]
        alpha = self.alpha.value() / 100
        show_over = alpha > 0 and self._over_lim is not None and self._over_lim[1] > self._over_lim[0]
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
                ax.set_title(title, fontsize=9)
                ax.set_xlim(-0.5, shape2d[0] - 0.5)
                ax.set_ylim(-0.5, shape2d[1] - 0.5)
                self._artists.append({"anat": anat_im, "over": over_im, "k": None, "marks": []})
        for ax, art, (title, plane, (h, w)) in zip(self.axes, self._artists, VIEWS):
            k = int(np.clip(round(c[plane]), 0, data.shape[plane] - 1))
            if art["k"] != k:
                sl = np.take(data, k, axis=plane)
                art["anat"].set_data(sl.T)
                art["k"] = k
                art["shape"] = sl.shape
            over = self._overlay_slice(plane, k, art["shape"]) if show_over else None
            if over is not None:
                art["over"].set_data(over.T)
                art["over"].set_cmap(self.cmap.currentText())
                art["over"].set_clim(*self._over_lim)
                art["over"].set_alpha(alpha)
                art["over"].set_visible(True)
            else:
                art["over"].set_visible(False)
            for line in art["marks"]:
                line.remove()
            art["marks"] = []
            for ras, marker, colour in self.markers:
                m = self.anat.ras_to_canon(ras)[0]
                near = abs(m[plane] - k) <= 2
                art["marks"] += ax.plot(m[h], m[w], marker, color=colour, markersize=8 if near else 6,
                                        markeredgewidth=1.5, alpha=1.0 if near else 0.35,
                                        linestyle="none")
        self.canvas.draw_idle()

    def _scrolled(self, event) -> None:
        if self.anat is None or event.inaxes not in self.axes:
            return
        _, plane, _ = VIEWS[self.axes.index(event.inaxes)]
        c = self.anat.ras_to_canon(self.cursor)[0]
        c[plane] += 1 if event.button == "up" else -1
        self.cursor = self.anat.canon_to_ras(c[None])[0]
        self._request_draw()
