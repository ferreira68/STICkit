#!/usr/bin/env python
"""Generate Mermaid call graph using PyCG when available, else a pure-AST fallback.
Writes:
  - docs/_mermaid/callgraph.json  (PyCG adjacency list when PyCG used)
  - docs/_mermaid/callgraph.mmd   (Mermaid graph TD)
"""

from __future__ import annotations

import ast
import json
import pathlib
import shutil
import subprocess  # nosec: B404
import sys
from subprocess import CompletedProcess
from typing import (
    Dict,
    Iterable,
    List,
    Sequence,
    Set,
    Tuple,
    cast,
)


REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "stickit"
OUT_DIR = REPO_ROOT / "docs" / "_mermaid"
JSON_OUT = OUT_DIR / "callgraph.json"
MMD_OUT = OUT_DIR / "callgraph.mmd"


def chunked(seq: List[str], size: int) -> Iterable[List[str]]:
    for i in range(0, len(seq), size):
        yield seq[i : i + size]


def list_py_files() -> list[str]:
    return [str(p) for p in SRC_DIR.rglob("*.py")]


def run(cmd: Sequence[str] | str) -> str:
    proc: CompletedProcess[str] = subprocess.run(cmd, capture_output=True, text=True, timeout=30)  # nosec B603
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr or proc.stdout)
    return proc.stdout


def merge_adj(a: Dict[str, list[str]], b: Dict[str, list[str]]) -> Dict[str, list[str]]:
    out = {k: list(v) for k, v in a.items()}
    for k, v in b.items():
        out.setdefault(k, [])
        seen = set(out[k])
        for x in v:
            if x not in seen:
                out[k].append(x)
                seen.add(x)
    return out


def try_pycg(files: list[str]) -> Dict[str, list[str]] | None:
    """Run PyCG on a possibly large set of files, chunking to avoid arg limits."""
    # Prefer module form to ensure Poetry venv is used
    if shutil.which("pycg") is None:
        return None
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    merged: Dict[str, list[str]] = {}
    chunk_size = 400
    for idx, chunk in enumerate(chunked(files, chunk_size), 1):
        tmp_json = OUT_DIR / f"callgraph.tmp.{idx}.json"
        cmd = [
            sys.executable,
            "-m",
            "pycg",
            "--package",
            "stickit",
            *chunk,
            "-o",
            str(tmp_json),
        ]
        try:
            run(cmd)
            if tmp_json.exists():
                part = json.loads(tmp_json.read_text(encoding="utf-8") or "{}")
                b_typed: Dict[str, list[str]] = cast(Dict[str, list[str]], part) if isinstance(part, dict) else {}
                merged = merge_adj(merged, b_typed)
                tmp_json.unlink(missing_ok=True)
        except Exception:
            return None
    return merged


# -------------------- AST fallback (dependency-free) -------------------------


def module_name(py_path: pathlib.Path) -> str:
    rel = py_path.relative_to(REPO_ROOT)
    return ".".join(rel.with_suffix("").parts)


def dotted_name(node: ast.AST) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        parts = []
        cur: ast.AST | None = node
        while isinstance(cur, ast.Attribute):
            parts.append(cur.attr)
            cur = cur.value
        if isinstance(cur, ast.Name):
            parts.append(cur.id)
            return ".".join(reversed(parts))
    return None


class CallGraphVisitor(ast.NodeVisitor):
    def __init__(self, mod: str) -> None:
        self.mod = mod
        self.stack: list[str] = []
        self.edges: Set[Tuple[str, str]] = set()

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:  # noqa: N802
        self.stack.append(node.name)
        self.generic_visit(node)
        self.stack.pop()

    # mypy: handle async separately to avoid incompatible arg type
    def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:  # noqa: N802
        self.stack.append(node.name)
        self.generic_visit(node)
        self.stack.pop()

    def visit_ClassDef(self, node: ast.ClassDef) -> None:  # noqa: N802
        self.stack.append(node.name)
        self.generic_visit(node)
        self.stack.pop()

    def visit_Call(self, node: ast.Call) -> None:  # noqa: N802
        src = self._qualname(self.stack[-1] if self.stack else "<module>")
        callee = dotted_name(node.func)
        if callee:
            self.edges.add((src, callee))
        self.generic_visit(node)

    def _qualname(self, leaf: str) -> str:
        if not self.stack:
            return f"{self.mod}.{leaf}"
        return f"{self.mod}." + ".".join(self.stack[:-1] + [leaf])


def ast_callgraph(files: list[str]) -> Dict[str, list[str]]:
    edges: Set[Tuple[str, str]] = set()
    for f in files:
        p = pathlib.Path(f)
        mod = module_name(p)
        try:
            src = p.read_text(encoding="utf-8")
            tree = ast.parse(src, filename=str(p))
        except Exception:
            continue
        v = CallGraphVisitor(mod)
        v.visit(tree)
        edges |= v.edges

    adj: Dict[str, list[str]] = {}
    for a, b in edges:
        adj.setdefault(a, []).append(b)
    return adj


# -------------------- Mermaid writer ----------------------------------------


def to_mermaid(adj: Dict[str, list[str]]) -> str:
    nodes: Set[str] = set()
    for a, dests in adj.items():
        nodes.add(a)
        nodes |= set(dests)
    lines = ["graph TD"]
    for n in sorted(nodes):
        safe = n.replace('"', "'")
        lines.append(f'  "{safe}"["{safe}"]')
    for a, dests in sorted(adj.items()):
        a_s = a.replace('"', "'")
        for b in sorted(set(dests)):
            b_s = b.replace('"', "'")
            lines.append(f'  "{a_s}" --> "{b_s}"')
    return "\n".join(lines)


# -------------------- main ---------------------------------------------------


def main() -> None:
    files = list_py_files()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    adj = try_pycg(files)
    if adj is None:
        adj = ast_callgraph(files)

    try:
        JSON_OUT.write_text(json.dumps(adj, indent=2), encoding="utf-8")
    except Exception as e:
        print(f"[callgraph] JSON write skipped: {e}", file=sys.stderr)

    MMD_OUT.write_text(to_mermaid(adj), encoding="utf-8")
    print(f"Wrote {MMD_OUT} (nodes={len(adj)})")


if __name__ == "__main__":
    main()
