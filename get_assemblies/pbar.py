#!/usr/bin/env python
# -*- coding: utf-8 -*-
from rich.console import Console
from rich.style import StyleType
from rich.progress import (
    TextColumn,
    SpinnerColumn,
    ProgressType,
    ProgressColumn,
    BarColumn,
    TimeRemainingColumn,
    Progress,
)
from typing import (
    Union,
    Sequence,
    Iterable,
    Optional,
    Callable,
    List,
)


def track_wide(
    sequence: Union[Sequence[ProgressType], Iterable[ProgressType]],
    description="Working...",
    total: int = None,
    auto_refresh=True,
    console: Optional[Console] = None,
    transient: bool = False,
    get_time: Callable[[], float] = None,
    refresh_per_second: float = None,
    style: StyleType = "bar.back",
    complete_style: StyleType = "bar.complete",
    finished_style: StyleType = "bar.finished",
    pulse_style: StyleType = "bar.pulse",
    update_period: float = 0.1,
    disable: bool = False,
) -> Iterable[ProgressType]:

    columns: List["ProgressColumn"]

    columns = [SpinnerColumn()]
    columns.extend(
        (
            TextColumn("[progress.description]{task.description}"),
            BarColumn(
                style=style,
                complete_style=complete_style,
                finished_style=finished_style,
                pulse_style=pulse_style,
                bar_width=None,
            ),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeRemainingColumn(),
        )
    )
    progress = Progress(
        *columns,
        auto_refresh=auto_refresh,
        console=console,
        transient=transient,
        get_time=get_time,
        refresh_per_second=refresh_per_second,
        disable=disable,
    )

    with progress:
        yield from progress.track(
            sequence,
            total=total,
            description=description.rjust(17, ' '),
            update_period=update_period,
        )
