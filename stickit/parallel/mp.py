from __future__ import annotations
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Iterable, Callable, Any


def run_parallel(
    items: Iterable[Any], fn: Callable[[Any], Any], max_workers: int | None = None
) -> list[Any]:
    """Run CPU-bound tasks in parallel using processes."""
    out: list[Any] = []
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = [ex.submit(fn, x) for x in items]
        for f in as_completed(futures):
            out.append(f.result())
    return out
