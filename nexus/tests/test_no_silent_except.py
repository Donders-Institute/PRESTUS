"""No silently swallowed exceptions in the Nexus without a stated reason.

An except clause is flagged when it is bare or catches Exception/BaseException, or when its
body only passes or continues. A flagged clause needs a marker comment on its except line:

    except OSError:  # best-effort: <why losing this error is harmless>
    except Exception as exc:  # reported: <where the user sees it>

Clauses that re-raise are not flagged: they do not swallow.
"""
from __future__ import annotations

import ast
import textwrap
from pathlib import Path

PACKAGE = Path(__file__).resolve().parent.parent / "prestus_nexus"
MARKERS = ("# best-effort:", "# reported:")
BROAD = {"Exception", "BaseException"}


def _is_broad(node: ast.ExceptHandler) -> bool:
    t = node.type
    if t is None:
        return True
    names = t.elts if isinstance(t, ast.Tuple) else [t]
    return any(isinstance(n, ast.Name) and n.id in BROAD for n in names)


def _is_silent(node: ast.ExceptHandler) -> bool:
    return all(isinstance(s, (ast.Pass, ast.Continue))
               or (isinstance(s, ast.Expr) and isinstance(s.value, ast.Constant))
               for s in node.body)


def _reraises(node: ast.ExceptHandler) -> bool:
    return any(isinstance(s, ast.Raise) for s in node.body)


def violations(source: str, name: str = "<src>") -> list[str]:
    lines = source.splitlines()
    out = []
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, ast.ExceptHandler) or _reraises(node):
            continue
        if not (_is_broad(node) or _is_silent(node)):
            continue
        line = lines[node.lineno - 1]
        marker = next((m for m in MARKERS if m in line), None)
        if marker is None or not line.split(marker, 1)[1].strip():
            out.append(f"{name}:{node.lineno}: {line.strip()}")
    return out


def test_every_broad_or_silent_except_states_why():
    found = []
    files = sorted(PACKAGE.rglob("*.py"))
    assert files, f"no sources under {PACKAGE}"
    for path in files:
        found += violations(path.read_text(encoding="utf-8"), str(path.relative_to(PACKAGE)))
    assert not found, ("except clauses that swallow errors without a reason; narrow them, "
                       "report the error, or add '# best-effort: <reason>' / '# reported: <where>' "
                       "on the except line:\n" + "\n".join(found))


def test_checker_flags_each_kind():          # negative control: the scan does find offenders
    src = textwrap.dedent("""
        try: f()
        except: g()
        try: f()
        except Exception: g()
        try: f()
        except (ValueError, BaseException): g()
        try: f()
        except OSError: pass
        for x in y:
            try: f()
            except KeyError: continue
        try: f()
        except OSError:  # best-effort:
            pass
    """)
    assert len(violations(src)) == 6


def test_checker_accepts_reasons_narrow_handling_and_reraise():
    src = textwrap.dedent("""
        try: f()
        except OSError:  # best-effort: temp file cleanup
            pass
        try: f()
        except Exception as exc:  # reported: in the status line
            show(exc)
        try: f()
        except ValueError:
            log()
        try: f()
        except BaseException:
            cleanup()
            raise
    """)
    assert violations(src) == []
