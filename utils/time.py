# これと同じです: https://github.com/gnya/toon/blob/main/toon/utils/time.py
from __future__ import annotations

import inspect
import time
from inspect import Parameter
from types import TracebackType
from typing import Any, Callable


class Stopwatch:
    def __init__(
        self, name: str, logger: Callable[[str], None] = print, count: int = 0
    ):
        self.name = name
        self.logger = logger
        self.count = count

    def start(self):
        self.time_start = time.perf_counter()

    def stop(self):
        self.time_end = time.perf_counter()
        mils = (self.time_end - self.time_start) * 1000

        self.logger(f"{self.name} {self.count}: {round(mils, 1)} ms")

    def lap(self) -> Stopwatch:
        self.stop()
        instance = Stopwatch(self.name, self.logger, self.count + 1)
        instance.start()

        return instance

    def __enter__(self):
        self.start()

    def __exit__(
        self,
        type: type[BaseException] | None,
        value: BaseException | None,
        traceback: TracebackType | None,
    ):
        self.stop()


def time_start(name: str, logger: Callable[[str], None] = print) -> Stopwatch:
    timer = Stopwatch(name, logger)
    timer.start()

    return timer


def timeit(
    args: Any = None, *, loops: int = 1, logger: Callable[[str], None] = print
) -> Callable[..., Any]:
    def wrapper(func: Callable[..., Any]) -> Callable[..., Any]:
        signature = inspect.signature(func)

        namespace = {"Stopwatch": Stopwatch, "func": func, "logger": logger}
        params = []

        for k, p in signature.parameters.items():
            if p.default is Parameter.empty:
                params.append(k)
            else:
                namespace[f"{k}_default"] = p.default
                params.append(f"{k} = {k}_default")

        params = ", ".join(params)
        assigns = ", ".join(signature.parameters.keys())

        code = (
            f"def wrapper({params}):\n"
            + "    with Stopwatch(func.__name__, logger):\n"
            + f"        if {loops > 1}:\n"
            + f"            for _ in range({loops - 1}):\n"
            + f"                func({assigns})\n"
            + "\n"
            + f"        return func({assigns})"
        )

        exec(code, namespace)

        return namespace["wrapper"]

    return wrapper(args) if callable(args) else wrapper
