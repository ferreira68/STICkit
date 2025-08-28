from __future__ import annotations
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Iterable, Callable, Any

def run_parallel(items: Iterable[Any], fn: Callable[[Any], Any], max_workers: int | None = None) -> list[Any]:
    """Run CPU-bound tasks in parallel using processes.

    Args:
        items: Iterable of picklable inputs.
        fn: Top-level picklable function.
        max_workers: Number of processes. Defaults to os.cpu_count().

    Returns:
        List of results in completion order.
    """
    out: list[Any] = []
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = [ex.submit(fn, x) for x in items]
        for f in as_completed(futures):
            out.append(f.result())
    return out
