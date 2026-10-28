"""Study manifest: one row per run (subject [+ condition]) with per-run config columns.

Same format as PRESTUS' study_read_manifest.m: tab- (.tsv) or comma-separated
(.csv), reserved columns subject_id / condition / include / config / notes, and
any other column a dotted config key ('placement.localite.target',
'transducer(2).serial'). Empty cells keep the study value.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

RESERVED = ("subject_id", "condition", "include", "config", "notes")
_IGNORED = {"notes", "note", "comment"}
_KEY_RE = re.compile(r"^[A-Za-z]\w*(\(\d+\))?(\.[A-Za-z]\w*(\(\d+\))?)*$")
_SUBJ_RE = re.compile(r"^(?:sub-)?(\d+)$")
_TRUE = {"", "1", "true", "yes", "y"}
_FALSE = {"0", "false", "no", "n"}


@dataclass
class Issue:
    level: str          # 'error' | 'warning'
    row: int | None     # 0-based data row, None = whole manifest
    column: str | None
    message: str

    def where(self) -> str:
        parts = []
        if self.row is not None:
            parts.append(f"row {self.row + 1}")
        if self.column:
            parts.append(self.column)
        return ", ".join(parts) or "manifest"


@dataclass
class Manifest:
    columns: list[str] = field(default_factory=lambda: ["subject_id", "condition"])
    rows: list[dict[str, str]] = field(default_factory=list)
    path: Path | None = None

    # ── io ──────────────────────────────────────────────────────────────
    @classmethod
    def load(cls, path: str | Path) -> "Manifest":
        path = Path(path)
        delim = "\t" if path.suffix.lower() == ".tsv" else ","
        lines = [ln for ln in path.read_text(encoding="utf-8").splitlines()
                 if ln.strip() and not ln.lstrip().startswith("#")]
        reader = csv.reader(lines, delimiter=delim)
        table = [list(r) for r in reader]
        if not table:
            return cls(path=path)
        columns = [c.strip() for c in table[0]]
        rows = []
        for cells in table[1:]:
            cells = [c.strip() for c in cells] + [""] * (len(columns) - len(cells))
            rows.append(dict(zip(columns, cells[:len(columns)])))
        return cls(columns=columns, rows=rows, path=path)

    def save(self, path: str | Path | None = None) -> Path:
        path = Path(path or self.path)
        delim = "\t" if path.suffix.lower() == ".tsv" else ","
        with path.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter=delim, lineterminator="\n")
            w.writerow(self.columns)
            for row in self.rows:
                w.writerow([row.get(c, "") for c in self.columns])
        self.path = path
        return path

    # ── editing ─────────────────────────────────────────────────────────
    def key_columns(self) -> list[str]:
        return [c for c in self.columns if c.lower() not in RESERVED and c.lower() not in _IGNORED]

    def add_column(self, name: str) -> None:
        if name not in self.columns:
            self.columns.append(name)

    def remove_column(self, name: str) -> None:
        if name == "subject_id":
            raise ValueError("subject_id is required")
        self.columns = [c for c in self.columns if c != name]
        for row in self.rows:
            row.pop(name, None)

    def add_row(self, values: dict[str, str] | None = None) -> None:
        values = dict(values or {})
        for k in values:
            self.add_column(k)
        self.rows.append({c: str(values.get(c, "")) for c in self.columns})

    # ── derived ─────────────────────────────────────────────────────────
    def affixes(self, study_affix: str = "") -> list[str]:
        """Output affix per row, as study_run_table.m derives it."""
        base = re.sub(r"^_desc-", "", study_affix or "")
        out = []
        for row in self.rows:
            explicit = row.get("io.output_affix", "")
            cond = row.get("condition", "")
            label = explicit or (f"{base}-{cond}" if base and cond else cond or base)
            out.append(sanitize_affix(label))
        return out

    def issues(self, schema=None, study_affix: str = "") -> list[Issue]:
        """Problems study_read_manifest / study_run_table would report, plus unknown keys."""
        issues: list[Issue] = []
        if "subject_id" not in [c.lower() for c in self.columns]:
            issues.append(Issue("error", None, None, "the manifest needs a subject_id column"))
        for c in self.key_columns():
            if not _KEY_RE.match(c):
                issues.append(Issue("error", None, c, "not a reserved column nor a dotted config key"))
            elif schema is not None and re.sub(r"\(\d+\)", "", c) not in schema:
                hint = schema.suggest(re.sub(r"\(\d+\)", "", c))
                issues.append(Issue("warning", None, c, "not a key of config_default.yaml"
                                    + (f" — did you mean {hint}?" if hint else "")))
        seen: dict[tuple[int, str], int] = {}
        affixes = self.affixes(study_affix)
        for i, row in enumerate(self.rows):
            m = _SUBJ_RE.match(row.get("subject_id", "").strip())
            if not m:
                issues.append(Issue("error", i, "subject_id", "must be an integer or sub-NNN"))
                continue
            inc = row.get("include", "").strip().lower()
            if inc not in _TRUE | _FALSE:
                issues.append(Issue("error", i, "include", "use 1/0, true/false or yes/no"))
            if inc in _FALSE:
                continue
            key = (int(m.group(1)), affixes[i])
            if key in seen:
                issues.append(Issue("error", i, "condition",
                                    f"same outputs as row {seen[key] + 1}; give the rows distinct conditions"))
            else:
                seen[key] = i
        return issues

    def conditions(self) -> list[str]:
        return sorted({r.get("condition", "") for r in self.rows})

    def subjects(self) -> list[int]:
        ids = set()
        for r in self.rows:
            m = _SUBJ_RE.match(r.get("subject_id", "").strip())
            if m:
                ids.add(int(m.group(1)))
        return sorted(ids)


def parse_subject_ids(text: str) -> list[str]:
    """Subject IDs from '1-12, 15 sub-020' (ranges, commas or spaces); zero padding kept.

    Raises ValueError naming the first token that is not an ID or a range."""
    ids: list[str] = []
    for token in re.split(r"[,;\s]+", text.strip()):
        if not token:
            continue
        rng = re.fullmatch(r"(?:sub-)?(\d+)-(?:sub-)?(\d+)", token)
        if rng:
            a, b = rng.groups()
            if int(b) < int(a):
                raise ValueError(f"'{token}': the range runs backwards")
            width = len(a) if a.startswith("0") else 0
            ids += [str(i).zfill(width) for i in range(int(a), int(b) + 1)]
        elif _SUBJ_RE.match(token):
            ids.append(token)
        else:
            raise ValueError(f"'{token}' is not a subject ID (use 3, 003, sub-003 or 1-12)")
    return list(dict.fromkeys(ids))


def sanitize_affix(label: str) -> str:
    """Python twin of sanitize_output_affix.m: '_desc-<label>', a leading 'desc-' not repeated."""
    if not label or not isinstance(label, str):
        return ""
    affix = re.sub(r"[^a-zA-Z0-9_-]", "_", label)
    if affix.startswith("_desc-"):
        return affix
    affix = re.sub(r"^desc-", "", affix.lstrip("_"))
    return "_desc-" + affix if affix else ""


def resolve_manifest(config_path: Path | None, value: str | None) -> Path | None:
    """study.manifest resolved against the study config's folder."""
    if not value:
        return None
    p = Path(value).expanduser()
    if not p.is_absolute() and config_path is not None:
        p = Path(config_path).parent / p
    return p


def from_neuronav_summary(summary_csv: str | Path, targets: list[str] | None = None) -> Manifest:
    """One localite run per ingested (subject, session, target) of neuronav_ingest_summary.csv.

    condition = target, suffixed with the session when a subject has the same
    target in several sessions. See neuronav_import for per-session
    (multi-transducer) grouping and the JSON sidecars.
    """
    from .neuronav_import import Entry, build_manifest
    with Path(summary_csv).open(encoding="utf-8", newline="") as f:
        rows = [r for r in csv.DictReader(f) if (r.get("target") or "").strip()]
    if targets:
        rows = [r for r in rows if r["target"] in targets]
    entries = sorted({Entry(r["sub"].strip(), r["ses"].strip(), r["target"].strip(), "summary") for r in rows},
                     key=lambda e: (e.sub, e.ses, e.target))
    return build_manifest(entries, "target")[0]
