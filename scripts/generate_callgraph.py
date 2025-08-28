#!/usr/bin/env python
"""Generate a Mermaid call graph from pyan3 output (non-breaking).
Writes docs/_mermaid/callgraph.mmd
"""
from __future__ import annotations
import subprocess, re, pathlib, sys

repo_root = pathlib.Path(__file__).resolve().parents[1]
src_dir = repo_root / "stickit"
out_file = repo_root / "docs" / "_mermaid" / "callgraph.mmd"

def run(cmd):
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr, file=sys.stderr)
        sys.exit(proc.returncode)
    return proc.stdout

try:
    dot = run([sys.executable, "-m", "pyan3", str(src_dir), "--no-defines", "--colored", "--grouped", "--dot"])
except FileNotFoundError:
    print("pyan3 is not installed. `poetry add -G dev pyan3` or merge supplement.")
    sys.exit(1)

edges = []
for line in dot.splitlines():
    m = re.search(r'\"([^\"]+)\"\s*->\s*\"([^\"]+)\"', line)
    if m:
        a, b = m.group(1), m.group(2)
        a_s = a.split(".")[-1]
        b_s = b.split(".")[-1]
        edges.append((a_s, b_s))

nodes = sorted({n for e in edges for n in e})
mermaid = ["graph TD"]
for n in nodes:
    mermaid.append(f"  {n}[{n}]")
for a, b in edges:
    mermaid.append(f"  {a} --> {b}")

out_file.parent.mkdir(parents=True, exist_ok=True)
out_file.write_text("\n".join(mermaid))
print(f"Wrote {out_file}")
