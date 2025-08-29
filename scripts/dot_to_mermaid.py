#!/usr/bin/env python
import re
import sys
import pathlib

src = pathlib.Path(sys.argv[1])
dst = pathlib.Path(sys.argv[2])
text = src.read_text(encoding="utf-8")
edges = re.findall(r'"([^"]+)"\s*->\s*"([^"]+)"', text)
nodes = sorted({n for e in edges for n in e})
lines = ["graph TD"]
for n in nodes:
    safe = n.replace('"', "'")
    lines.append(f'  "{safe}"["{safe}"]')
for a, b in edges:
    a = a.replace('"', "'")
    b = b.replace('"', "'")
    lines.append(f'  "{a}" --> "{b}"')
dst.write_text("\n".join(lines), encoding="utf-8")
