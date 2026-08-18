"""Actuarial reserve modeling with reusable ODE solvers.

This script validates Euler, Heun, and RK4 methods against closed-form ODE
solutions, then applies the same solver interface to Thiele's differential
equation for a simplified continuous life insurance reserve model.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from html import escape
from pathlib import Path
from typing import Callable, Iterable


State = float
Derivative = Callable[[float, State], State]
Solver = Callable[[Derivative, State, Iterable[float]], list[State]]


def euler(f: Derivative, y0: State, grid: Iterable[float]) -> list[State]:
    times = list(grid)
    values = [y0]
    for current, nxt in zip(times, times[1:]):
        h = nxt - current
        values.append(values[-1] + h * f(current, values[-1]))
    return values


def heun(f: Derivative, y0: State, grid: Iterable[float]) -> list[State]:
    times = list(grid)
    values = [y0]
    for current, nxt in zip(times, times[1:]):
        h = nxt - current
        predictor = values[-1] + h * f(current, values[-1])
        corrected_slope = 0.5 * (f(current, values[-1]) + f(nxt, predictor))
        values.append(values[-1] + h * corrected_slope)
    return values


def rk4(f: Derivative, y0: State, grid: Iterable[float]) -> list[State]:
    times = list(grid)
    values = [y0]
    for current, nxt in zip(times, times[1:]):
        h = nxt - current
        y = values[-1]
        k1 = f(current, y)
        k2 = f(current + 0.5 * h, y + 0.5 * h * k1)
        k3 = f(current + 0.5 * h, y + 0.5 * h * k2)
        k4 = f(nxt, y + h * k3)
        values.append(y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4))
    return values


def make_grid(start: float, end: float, steps: int) -> list[float]:
    if steps <= 0:
        raise ValueError("steps must be positive")
    h = (end - start) / steps
    return [start + i * h for i in range(steps + 1)]


@dataclass(frozen=True)
class ExponentialDecayModel:
    rate: float

    def derivative(self, _: float, y: float) -> float:
        return -self.rate * y

    def exact(self, t: float, y0: float) -> float:
        return y0 * math.exp(-self.rate * t)


@dataclass(frozen=True)
class LogisticGrowthModel:
    growth_rate: float
    carrying_capacity: float

    def derivative(self, _: float, y: float) -> float:
        return self.growth_rate * y * (1 - y / self.carrying_capacity)

    def exact(self, t: float, y0: float) -> float:
        ratio = (self.carrying_capacity - y0) / y0
        return self.carrying_capacity / (1 + ratio * math.exp(-self.growth_rate * t))


@dataclass(frozen=True)
class LifeReserveParams:
    interest_rate: float = 0.035
    mortality_rate: float = 0.012
    annual_premium: float = 1_600.0
    death_benefit: float = 100_000.0


@dataclass(frozen=True)
class ThieleReserveModel:
    params: LifeReserveParams

    def derivative(self, _: float, reserve: float) -> float:
        p = self.params
        return p.interest_rate * reserve + p.annual_premium - p.mortality_rate * (
            p.death_benefit - reserve
        )


SOLVERS: dict[str, Solver] = {
    "Euler": euler,
    "Heun": heun,
    "RK4": rk4,
}


def max_absolute_error(values: list[float], expected: list[float]) -> float:
    return max(abs(a - b) for a, b in zip(values, expected))


def validate_solvers() -> list[dict[str, float | str]]:
    cases = [
        ("Exponential decay", ExponentialDecayModel(rate=0.35), 100.0, 0.0, 12.0, 120),
        ("Logistic growth", LogisticGrowthModel(0.45, 1_000.0), 80.0, 0.0, 20.0, 200),
    ]
    rows: list[dict[str, float | str]] = []
    for name, model, y0, start, end, steps in cases:
        grid = make_grid(start, end, steps)
        expected = [model.exact(t, y0) for t in grid]
        for solver_name, solver in SOLVERS.items():
            values = solver(model.derivative, y0, grid)
            rows.append(
                {
                    "model": name,
                    "solver": solver_name,
                    "steps": steps,
                    "max_absolute_error": max_absolute_error(values, expected),
                }
            )
    return rows


def run_reserve_scenarios() -> tuple[list[float], dict[str, list[float]]]:
    scenarios = {
        "Base": LifeReserveParams(),
        "Higher mortality": LifeReserveParams(mortality_rate=0.018),
        "Lower interest": LifeReserveParams(interest_rate=0.02),
        "Higher benefit": LifeReserveParams(death_benefit=125_000.0),
    }
    grid = make_grid(0.0, 30.0, 360)
    reserve_paths = {
        label: rk4(ThieleReserveModel(params).derivative, 0.0, grid)
        for label, params in scenarios.items()
    }
    return grid, reserve_paths


def write_solver_validation(rows: list[dict[str, float | str]], output_dir: Path) -> None:
    path = output_dir / "solver_validation_summary.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["model", "solver", "steps", "max_absolute_error"]
        )
        writer.writeheader()
        writer.writerows(rows)


def write_reserve_scenarios(
    grid: list[float], reserve_paths: dict[str, list[float]], output_dir: Path
) -> None:
    path = output_dir / "life_reserve_scenarios.csv"
    fieldnames = ["year", *reserve_paths.keys()]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for idx, year in enumerate(grid):
            row = {"year": year}
            row.update({label: values[idx] for label, values in reserve_paths.items()})
            writer.writerow(row)


def write_solver_validation_svg(rows: list[dict[str, float | str]], figure_dir: Path) -> None:
    labels = [f"{row['model']}\n{row['solver']}" for row in rows]
    errors = [float(row["max_absolute_error"]) for row in rows]
    colors = ["#3B6EA8", "#E2A53A", "#8B9A46"] * 2

    width, height = 980, 560
    left, top, right, bottom = 82, 74, 36, 132
    plot_w = width - left - right
    plot_h = height - top - bottom
    min_log = math.floor(math.log10(min(errors)))
    max_log = math.ceil(math.log10(max(errors)))

    def y_pos(value: float) -> float:
        scaled = (math.log10(value) - min_log) / (max_log - min_log)
        return top + plot_h * (1 - scaled)

    bar_gap = 18
    bar_w = (plot_w - bar_gap * (len(rows) - 1)) / len(rows)
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#FFFFFF"/>',
        '<text x="40" y="38" font-family="Arial, sans-serif" font-size="24" font-weight="700" fill="#222222">Solver Validation Against Closed-Form ODE Solutions</text>',
        '<text x="40" y="62" font-family="Arial, sans-serif" font-size="13" fill="#555555">Maximum absolute error on a logarithmic scale</text>',
    ]

    for tick in range(min_log, max_log + 1):
        y = y_pos(10**tick)
        svg.append(f'<line x1="{left}" y1="{y:.2f}" x2="{width-right}" y2="{y:.2f}" stroke="#E1E1E1"/>')
        svg.append(f'<text x="{left-10}" y="{y+4:.2f}" text-anchor="end" font-family="Arial, sans-serif" font-size="12" fill="#555555">1e{tick}</text>')

    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#333333"/>')
    svg.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{width-right}" y2="{top+plot_h}" stroke="#333333"/>')

    for idx, (label, error, color) in enumerate(zip(labels, errors, colors)):
        x = left + idx * (bar_w + bar_gap)
        y = y_pos(error)
        bar_h = top + plot_h - y
        svg.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_w:.2f}" height="{bar_h:.2f}" fill="{color}" stroke="#2B2B2B" stroke-width="0.8"/>')
        svg.append(f'<text x="{x + bar_w / 2:.2f}" y="{y - 7:.2f}" text-anchor="middle" font-family="Arial, sans-serif" font-size="11" fill="#333333">{error:.1e}</text>')
        first, second = label.split("\n")
        svg.append(f'<text x="{x + bar_w / 2:.2f}" y="{top + plot_h + 24}" text-anchor="middle" font-family="Arial, sans-serif" font-size="11" fill="#333333">{escape(first)}</text>')
        svg.append(f'<text x="{x + bar_w / 2:.2f}" y="{top + plot_h + 40}" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" font-weight="700" fill="#333333">{escape(second)}</text>')

    svg.append("</svg>")
    (figure_dir / "solver_validation.svg").write_text("\n".join(svg))


def write_reserve_scenarios_svg(
    grid: list[float], reserve_paths: dict[str, list[float]], figure_dir: Path
) -> None:
    palette = {
        "Base": "#3B6EA8",
        "Higher mortality": "#CC6B49",
        "Lower interest": "#8B9A46",
        "Higher benefit": "#A45A7A",
    }

    width, height = 980, 590
    left, top, right, bottom = 82, 82, 210, 70
    plot_w = width - left - right
    plot_h = height - top - bottom
    values = [value for path in reserve_paths.values() for value in path]
    min_y = math.floor(min(values) / 1000) * 1000
    max_y = math.ceil(max(values) / 1000) * 1000
    min_x, max_x = min(grid), max(grid)

    def x_pos(year: float) -> float:
        return left + plot_w * ((year - min_x) / (max_x - min_x))

    def y_pos(value: float) -> float:
        return top + plot_h * (1 - (value - min_y) / (max_y - min_y))

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#FFFFFF"/>',
        '<text x="40" y="40" font-family="Arial, sans-serif" font-size="25" font-weight="700" fill="#222222">Life Insurance Reserve Scenarios</text>',
        '<text x="40" y="64" font-family="Arial, sans-serif" font-size="13" fill="#555555">30-year continuous reserve paths under mortality, interest, and benefit sensitivity assumptions</text>',
    ]

    for tick in range(0, 31, 5):
        x = x_pos(tick)
        svg.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top+plot_h}" stroke="#EFEFEF"/>')
        svg.append(f'<text x="{x:.2f}" y="{top+plot_h+26}" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#555555">{tick}</text>')

    step = max(1000, int((max_y - min_y) / 5))
    for tick in range(int(min_y), int(max_y) + step, step):
        y = y_pos(tick)
        svg.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}" stroke="#E1E1E1"/>')
        svg.append(f'<text x="{left-10}" y="{y+4:.2f}" text-anchor="end" font-family="Arial, sans-serif" font-size="12" fill="#555555">{tick:,.0f}</text>')

    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#333333"/>')
    svg.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#333333"/>')
    svg.append(f'<text x="{left + plot_w / 2}" y="{height-18}" text-anchor="middle" font-family="Arial, sans-serif" font-size="13" fill="#333333">Policy year</text>')
    svg.append(f'<text x="22" y="{top + plot_h / 2}" transform="rotate(-90 22 {top + plot_h / 2})" text-anchor="middle" font-family="Arial, sans-serif" font-size="13" fill="#333333">Reserve</text>')

    sampled_grid = grid[::20]
    for label, path in reserve_paths.items():
        sampled_path = path[::20]
        if sampled_grid[-1] != grid[-1]:
            sampled_grid = [*sampled_grid, grid[-1]]
            sampled_path = [*sampled_path, path[-1]]
        points = " ".join(
            f"{x_pos(t):.2f},{y_pos(v):.2f}" for t, v in zip(sampled_grid, sampled_path)
        )
        svg.append(f'<polyline points="{points}" fill="none" stroke="{palette[label]}" stroke-width="3" stroke-linejoin="round" stroke-linecap="round"/>')

    legend_x = left + plot_w + 38
    legend_y = top + 16
    for idx, label in enumerate(reserve_paths):
        y = legend_y + idx * 30
        svg.append(f'<line x1="{legend_x}" y1="{y}" x2="{legend_x+28}" y2="{y}" stroke="{palette[label]}" stroke-width="4" stroke-linecap="round"/>')
        svg.append(f'<text x="{legend_x+38}" y="{y+4}" font-family="Arial, sans-serif" font-size="13" fill="#333333">{escape(label)}</text>')

    svg.append("</svg>")
    (figure_dir / "life_reserve_scenarios.svg").write_text("\n".join(svg))


def main() -> None:
    root = Path(__file__).resolve().parent
    output_dir = root / "outputs"
    figure_dir = root / "figures"
    output_dir.mkdir(exist_ok=True)
    figure_dir.mkdir(exist_ok=True)

    validation_rows = validate_solvers()
    grid, reserve_paths = run_reserve_scenarios()

    write_solver_validation(validation_rows, output_dir)
    write_reserve_scenarios(grid, reserve_paths, output_dir)
    write_solver_validation_svg(validation_rows, figure_dir)
    write_reserve_scenarios_svg(grid, reserve_paths, figure_dir)

    print("Solver validation")
    for row in validation_rows:
        print(
            f"- {row['model']} | {row['solver']}: "
            f"max error = {float(row['max_absolute_error']):.4e}"
        )
    print("\nGenerated outputs:")
    print(f"- {output_dir / 'solver_validation_summary.csv'}")
    print(f"- {output_dir / 'life_reserve_scenarios.csv'}")
    print(f"- {figure_dir / 'solver_validation.svg'}")
    print(f"- {figure_dir / 'life_reserve_scenarios.svg'}")


if __name__ == "__main__":
    main()
