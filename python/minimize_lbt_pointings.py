#!/usr/bin/env python3
"""
Minimize the number of LBT pointings needed to cover a galaxy catalog.

The field is treated as a fixed-orientation square on the local tangent plane:
25 x 25 arcmin by default. The code generates the finite set of candidate
pointings that can be optimal, then solves the resulting minimum set-cover
problem exactly with branch and bound. A final polishing step recenters each
selected field, keeping its required galaxies covered while trying to include
additional galaxies that fit in the same field. If --max-pointings is supplied,
the code can switch to a capped greedy selection. The capped modes either anchor
the next field on the highest-priority target not already covered, or maximize
the number of targets above a priority threshold covered by each field.

Input catalog:
  - CSV file
  - RA and Dec columns in decimal degrees
  - radius column, in arcmin by default
  - optional priority column, such as log stellar mass
  - optional object ID/name column

Example:
  python minimize_lbt_pointings.py galaxies.csv \
      --ra-col ra_deg --dec-col dec_deg --radius-col radius --id-col name \
      --pointings-out lbt_pointings.csv \
      --assignments-out lbt_assignments.csv

This script uses only the Python standard library.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from dataclasses import dataclass
from typing import Iterable

ARC_MIN_PER_RAD = 180.0 * 60.0 / math.pi
RAD_PER_ARC_MIN = math.pi / (180.0 * 60.0)


def popcount(value: int) -> int:
    if hasattr(int, "bit_count"):
        return value.bit_count()
    return bin(value).count("1")


@dataclass(frozen=True)
class Target:
    index: int
    target_id: str
    ra_deg: float
    dec_deg: float
    radius_arcmin: float
    priority: float
    vector: tuple[float, float, float]


@dataclass(frozen=True)
class Basis:
    origin: tuple[float, float, float]
    east: tuple[float, float, float]
    north: tuple[float, float, float]


@dataclass(frozen=True)
class Candidate:
    ra_deg: float
    dec_deg: float
    vector: tuple[float, float, float]
    covers: frozenset[int]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Optimize fixed-size LBT pointings for a galaxy catalog."
    )
    parser.add_argument("catalog", help="Input CSV catalog.")
    parser.add_argument("--ra-col", default="ra", help="RA column in degrees.")
    parser.add_argument("--dec-col", default="dec", help="Dec column in degrees.")
    parser.add_argument(
        "--radius-col",
        default="radius",
        help="Galaxy radius column. Default: radius.",
    )
    parser.add_argument(
        "--radius-unit",
        choices=["arcmin", "arcsec", "deg"],
        default="arcmin",
        help="Unit for --radius-col. Default: arcmin.",
    )
    parser.add_argument(
        "--id-col",
        default=None,
        help="Optional object ID/name column. If omitted, row numbers are used.",
    )
    parser.add_argument(
        "--priority-col",
        default=None,
        help="Optional priority column, e.g. logMstar. Used with --max-pointings.",
    )
    parser.add_argument(
        "--priority-order",
        choices=["desc", "asc"],
        default="desc",
        help="Sort priority descending or ascending. Default: desc.",
    )
    parser.add_argument(
        "--selection-mode",
        choices=["anchor-priority", "max-count-above-threshold"],
        default="anchor-priority",
        help=(
            "Capped selection behavior. anchor-priority starts with the highest-priority "
            "uncovered target. max-count-above-threshold maximizes the number of "
            "uncovered targets above --priority-threshold in each field."
        ),
    )
    parser.add_argument(
        "--priority-threshold",
        type=float,
        default=None,
        help=(
            "Priority threshold for --selection-mode max-count-above-threshold, "
            "e.g. 8.5 for logMstar > 8.5."
        ),
    )
    parser.add_argument(
        "--max-pointings",
        type=int,
        default=None,
        help="Maximum number of pointings. If set, select fields by priority instead of covering every target.",
    )
    parser.add_argument(
        "--fov-arcmin",
        type=float,
        default=25.0,
        help="Square field side length in arcmin. Default: 25.",
    )
    parser.add_argument(
        "--pointings-out",
        default="lbt_pointings.csv",
        help="Output CSV for optimized pointing centers.",
    )
    parser.add_argument(
        "--assignments-out",
        default="lbt_assignments.csv",
        help="Output CSV assigning each target to one selected pointing.",
    )
    parser.add_argument(
        "--time-limit",
        type=float,
        default=None,
        help="Optional exact solver time limit in seconds.",
    )
    return parser.parse_args()


def radius_to_arcmin(value: float, unit: str) -> float:
    if unit == "arcmin":
        return value
    if unit == "arcsec":
        return value / 60.0
    if unit == "deg":
        return value * 60.0
    raise ValueError(f"Unsupported radius unit: {unit}")


def dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def norm(v: tuple[float, float, float]) -> float:
    return math.sqrt(dot(v, v))


def normalize(v: tuple[float, float, float]) -> tuple[float, float, float]:
    length = norm(v)
    if length == 0.0:
        raise ValueError("Cannot normalize a zero-length vector.")
    return (v[0] / length, v[1] / length, v[2] / length)


def vector_from_radec(ra_deg: float, dec_deg: float) -> tuple[float, float, float]:
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    cos_dec = math.cos(dec)
    return (
        cos_dec * math.cos(ra),
        cos_dec * math.sin(ra),
        math.sin(dec),
    )


def radec_from_vector(v: tuple[float, float, float]) -> tuple[float, float]:
    x, y, z = normalize(v)
    ra = math.degrees(math.atan2(y, x)) % 360.0
    dec = math.degrees(math.asin(max(-1.0, min(1.0, z))))
    return ra, dec


def basis_from_vector(origin: tuple[float, float, float]) -> Basis:
    ra_deg, dec_deg = radec_from_vector(origin)
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    east = (-math.sin(ra), math.cos(ra), 0.0)
    north = (
        -math.sin(dec) * math.cos(ra),
        -math.sin(dec) * math.sin(ra),
        math.cos(dec),
    )
    return Basis(origin=normalize(origin), east=normalize(east), north=normalize(north))


def angular_sep_arcmin(
    a: tuple[float, float, float], b: tuple[float, float, float]
) -> float:
    cos_sep = max(-1.0, min(1.0, dot(a, b)))
    return math.acos(cos_sep) * ARC_MIN_PER_RAD


def tangent_offsets_arcmin(
    vector: tuple[float, float, float], basis: Basis
) -> tuple[float, float]:
    """
    Gnomonic tangent-plane offsets, positive east and north, in arcmin.

    For a 25 arcmin field this is effectively exact for field placement, and it
    handles RA wrap naturally because all calculations are vector-based.
    """
    denom = dot(vector, basis.origin)
    if denom <= 0.0:
        raise ValueError("Target is more than 90 degrees from the tangent origin.")
    x = dot(vector, basis.east) / denom
    y = dot(vector, basis.north) / denom
    return x * ARC_MIN_PER_RAD, y * ARC_MIN_PER_RAD


def vector_from_tangent_offset(
    basis: Basis, x_arcmin: float, y_arcmin: float
) -> tuple[float, float, float]:
    x = x_arcmin * RAD_PER_ARC_MIN
    y = y_arcmin * RAD_PER_ARC_MIN
    vector = (
        basis.origin[0] + x * basis.east[0] + y * basis.north[0],
        basis.origin[1] + x * basis.east[1] + y * basis.north[1],
        basis.origin[2] + x * basis.east[2] + y * basis.north[2],
    )
    return normalize(vector)


def mean_origin(targets: list[Target], indices: Iterable[int]) -> tuple[float, float, float]:
    sx = sy = sz = 0.0
    first = None
    for i in indices:
        first = i if first is None else first
        v = targets[i].vector
        sx += v[0]
        sy += v[1]
        sz += v[2]
    if norm((sx, sy, sz)) < 1e-12:
        if first is None:
            raise ValueError("Cannot build an origin for an empty component.")
        return targets[first].vector
    return normalize((sx, sy, sz))


def load_catalog(
    path: str,
    ra_col: str,
    dec_col: str,
    radius_col: str,
    radius_unit: str,
    priority_col: str | None,
    id_col: str | None,
) -> list[Target]:
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("Input CSV has no header row.")
        missing = [
            col for col in [ra_col, dec_col, radius_col] if col not in reader.fieldnames
        ]
        if id_col is not None and id_col not in reader.fieldnames:
            missing.append(id_col)
        if priority_col is not None and priority_col not in reader.fieldnames:
            missing.append(priority_col)
        if missing:
            raise ValueError(
                f"Missing column(s): {missing}. Available columns: {reader.fieldnames}"
            )

        targets: list[Target] = []
        for row_index, row in enumerate(reader):
            ra_deg = float(row[ra_col])
            dec_deg = float(row[dec_col])
            radius_arcmin = radius_to_arcmin(float(row[radius_col]), radius_unit)
            if radius_arcmin < 0.0:
                raise ValueError(
                    f"Radius cannot be negative for row {row_index}: {radius_arcmin}"
                )
            if priority_col is None:
                priority = float(row_index)
            else:
                raw_priority = row[priority_col]
                priority = float(raw_priority) if raw_priority not in ("", "nan", "NaN") else -math.inf
            target_id = row[id_col] if id_col else str(row_index)
            targets.append(
                Target(
                    index=row_index,
                    target_id=target_id,
                    ra_deg=ra_deg,
                    dec_deg=dec_deg,
                    radius_arcmin=radius_arcmin,
                    priority=priority,
                    vector=vector_from_radec(ra_deg, dec_deg),
                )
            )

    if not targets:
        raise ValueError("Input catalog is empty.")
    return targets


def connected_components(targets: list[Target], fov_arcmin: float) -> list[list[int]]:
    """
    Split targets into independent groups.

    If two targets are farther apart than the field diagonal, they cannot share
    one pointing. Components of this graph can be optimized separately.
    """
    n_targets = len(targets)
    adjacency = [[] for _ in range(n_targets)]

    for i in range(n_targets):
        for j in range(i + 1, n_targets):
            max_same_field = (
                math.sqrt(2.0)
                * max(0.0, fov_arcmin - targets[i].radius_arcmin - targets[j].radius_arcmin)
            )
            sep = angular_sep_arcmin(targets[i].vector, targets[j].vector)
            if sep <= max_same_field + 1e-8:
                adjacency[i].append(j)
                adjacency[j].append(i)

    seen = [False] * n_targets
    components: list[list[int]] = []
    for start in range(n_targets):
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        component = []
        while stack:
            i = stack.pop()
            component.append(i)
            for j in adjacency[i]:
                if not seen[j]:
                    seen[j] = True
                    stack.append(j)
        components.append(sorted(component))

    return components


def covered_by_center(
    center_vector: tuple[float, float, float],
    targets: list[Target],
    indices: Iterable[int],
    half_size_arcmin: float,
) -> frozenset[int]:
    basis = basis_from_vector(center_vector)
    covered = []
    for i in indices:
        try:
            x, y = tangent_offsets_arcmin(targets[i].vector, basis)
        except ValueError:
            continue
        radius = targets[i].radius_arcmin
        if (
            abs(x) + radius <= half_size_arcmin + 1e-8
            and abs(y) + radius <= half_size_arcmin + 1e-8
        ):
            covered.append(i)
    return frozenset(covered)


def candidate_from_center(
    center_vector: tuple[float, float, float],
    targets: list[Target],
    indices: Iterable[int],
    half_size_arcmin: float,
) -> Candidate:
    covers = covered_by_center(center_vector, targets, indices, half_size_arcmin)
    ra_deg, dec_deg = radec_from_vector(center_vector)
    return Candidate(
        ra_deg=ra_deg,
        dec_deg=dec_deg,
        vector=center_vector,
        covers=covers,
    )


def generate_candidates_for_component(
    targets: list[Target],
    component: list[int],
    fov_arcmin: float,
) -> list[Candidate]:
    """
    Generate every candidate square that can matter for an optimum.

    In a local tangent plane, some optimum can be shifted until a left/right
    edge touches a target and a bottom/top edge touches a target. Those edge
    combinations define a finite candidate list.
    """
    half = fov_arcmin / 2.0
    origin = mean_origin(targets, component)
    basis = basis_from_vector(origin)

    offsets = {
        i: tangent_offsets_arcmin(targets[i].vector, basis)
        for i in component
    }
    left_edges = sorted(
        set(
            [offsets[i][0] - targets[i].radius_arcmin for i in component]
            + [
                offsets[i][0] + targets[i].radius_arcmin - fov_arcmin
                for i in component
            ]
        )
    )
    bottom_edges = sorted(
        set(
            [offsets[i][1] - targets[i].radius_arcmin for i in component]
            + [
                offsets[i][1] + targets[i].radius_arcmin - fov_arcmin
                for i in component
            ]
        )
    )

    candidates_by_cover: dict[tuple[int, ...], Candidate] = {}

    for i in component:
        if targets[i].radius_arcmin > half + 1e-8:
            raise ValueError(
                f"Target {targets[i].target_id!r} has radius "
                f"{targets[i].radius_arcmin:.6f} arcmin, larger than half the "
                f"{fov_arcmin:.6f} arcmin field. It cannot be fully covered."
            )
        center_vector = targets[i].vector
        covers = covered_by_center(center_vector, targets, component, half)
        if covers:
            ra_deg, dec_deg = radec_from_vector(center_vector)
            candidates_by_cover[tuple(sorted(covers))] = Candidate(
                ra_deg=ra_deg,
                dec_deg=dec_deg,
                vector=center_vector,
                covers=covers,
            )

    for left in left_edges:
        cx = left + half
        possible_x = {
            i
            for i in component
            if abs(offsets[i][0] - cx) + targets[i].radius_arcmin <= half + 1e-8
        }
        if not possible_x:
            continue
        for bottom in bottom_edges:
            cy = bottom + half
            possible_y = {
                i
                for i in component
                if abs(offsets[i][1] - cy) + targets[i].radius_arcmin <= half + 1e-8
            }
            if not possible_y or not (possible_x & possible_y):
                continue

            center_vector = vector_from_tangent_offset(basis, cx, cy)
            covers = covered_by_center(center_vector, targets, component, half)
            if not covers:
                continue
            ra_deg, dec_deg = radec_from_vector(center_vector)
            candidates_by_cover[tuple(sorted(covers))] = Candidate(
                ra_deg=ra_deg,
                dec_deg=dec_deg,
                vector=center_vector,
                covers=covers,
            )

    return list(candidates_by_cover.values())


def greedy_solution(universe: int, masks: list[int]) -> list[int]:
    uncovered = universe
    chosen = []
    while uncovered:
        best_j = max(
            range(len(masks)),
            key=lambda j: popcount(masks[j] & uncovered),
        )
        if (masks[best_j] & uncovered) == 0:
            raise RuntimeError("No candidate can cover the remaining targets.")
        chosen.append(best_j)
        uncovered &= ~masks[best_j]
    return chosen


def exact_set_cover(
    masks: list[int],
    n_targets: int,
    time_limit: float | None,
) -> list[int]:
    """Solve minimum set cover exactly with branch and bound."""
    start_time = time.monotonic()
    universe = (1 << n_targets) - 1
    max_cover = max(popcount(mask) for mask in masks)
    point_to_candidates = [[] for _ in range(n_targets)]

    for j, mask in enumerate(masks):
        for bit in range(n_targets):
            if mask & (1 << bit):
                point_to_candidates[bit].append(j)

    for bit, choices in enumerate(point_to_candidates):
        if not choices:
            raise RuntimeError(f"Internal error: local target {bit} has no candidates.")

    best = greedy_solution(universe, masks)
    memo: dict[int, int] = {}

    def timed_out() -> bool:
        return time_limit is not None and time.monotonic() - start_time > time_limit

    def search(uncovered: int, chosen: list[int]) -> None:
        nonlocal best
        if timed_out():
            raise TimeoutError(
                f"Exact solver exceeded {time_limit} seconds; "
                f"best solution so far used {len(best)} pointings."
            )
        if not uncovered:
            if len(chosen) < len(best):
                best = chosen[:]
            return
        if len(chosen) >= len(best):
            return

        lower_bound = math.ceil(popcount(uncovered) / max_cover)
        if len(chosen) + lower_bound >= len(best):
            return

        if memo.get(uncovered, 10**9) <= len(chosen):
            return
        memo[uncovered] = len(chosen)

        uncovered_bits = [
            bit for bit in range(n_targets) if uncovered & (1 << bit)
        ]
        pivot = min(
            uncovered_bits,
            key=lambda bit: sum(1 for j in point_to_candidates[bit] if masks[j] & uncovered),
        )
        options = [
            j for j in point_to_candidates[pivot] if masks[j] & uncovered
        ]
        options.sort(key=lambda j: popcount(masks[j] & uncovered), reverse=True)

        seen_next: set[int] = set()
        for j in options:
            next_uncovered = uncovered & ~masks[j]
            if next_uncovered in seen_next:
                continue
            seen_next.add(next_uncovered)
            chosen.append(j)
            search(next_uncovered, chosen)
            chosen.pop()

    search(universe, [])
    return best


def solve_component(
    component: list[int],
    candidates: list[Candidate],
    time_limit: float | None,
) -> list[int]:
    local_index = {global_i: local_i for local_i, global_i in enumerate(component)}
    mask_to_candidate: dict[int, int] = {}

    for j, candidate in enumerate(candidates):
        mask = 0
        for global_i in candidate.covers:
            mask |= 1 << local_index[global_i]
        if mask:
            mask_to_candidate.setdefault(mask, j)

    masks = list(mask_to_candidate.keys())
    candidate_indices = [mask_to_candidate[mask] for mask in masks]
    chosen_mask_indices = exact_set_cover(masks, len(component), time_limit)
    return [candidate_indices[j] for j in chosen_mask_indices]


def best_priority_candidate_for_anchor(
    anchor: int,
    targets: list[Target],
    covered: set[int],
    fov_arcmin: float,
) -> Candidate:
    """
    Choose the best field that includes the anchor target.

    The anchor is fixed by priority order. The field center is allowed to move
    anywhere that keeps the anchor covered, and is scored by the total priority
    of not-yet-covered targets that fit in the field.
    """
    half = fov_arcmin / 2.0
    if targets[anchor].radius_arcmin > half + 1e-8:
        raise ValueError(
            f"Target {targets[anchor].target_id!r} has radius "
            f"{targets[anchor].radius_arcmin:.6f} arcmin, larger than half the "
            f"{fov_arcmin:.6f} arcmin field. It cannot be fully covered."
        )

    anchor_vector = targets[anchor].vector
    basis = basis_from_vector(anchor_vector)
    local_indices = []
    offsets: dict[int, tuple[float, float]] = {}

    for i, target in enumerate(targets):
        if target.radius_arcmin > half + 1e-8:
            continue
        try:
            x, y = tangent_offsets_arcmin(target.vector, basis)
        except ValueError:
            continue
        max_relevant_offset = math.sqrt(2.0) * fov_arcmin + target.radius_arcmin
        if abs(x) <= max_relevant_offset and abs(y) <= max_relevant_offset:
            local_indices.append(i)
            offsets[i] = (x, y)

    left_edges = sorted(
        set(
            [offsets[i][0] - targets[i].radius_arcmin for i in local_indices]
            + [
                offsets[i][0] + targets[i].radius_arcmin - fov_arcmin
                for i in local_indices
            ]
        )
    )
    bottom_edges = sorted(
        set(
            [offsets[i][1] - targets[i].radius_arcmin for i in local_indices]
            + [
                offsets[i][1] + targets[i].radius_arcmin - fov_arcmin
                for i in local_indices
            ]
        )
    )

    all_indices = list(range(len(targets)))
    best = candidate_from_center(anchor_vector, targets, all_indices, half)

    def score(candidate: Candidate) -> tuple[float, int, int, float]:
        newly_covered = candidate.covers - covered
        priority_sum = sum(
            targets[i].priority for i in newly_covered if math.isfinite(targets[i].priority)
        )
        return (
            priority_sum,
            len(newly_covered),
            len(candidate.covers),
            min_edge_margin_arcmin(candidate.vector, targets, [anchor], half),
        )

    best_score = score(best)
    seen_centers: set[tuple[float, float]] = set()

    for left in left_edges:
        cx = left + half
        for bottom in bottom_edges:
            cy = bottom + half
            center_key = (round(cx, 8), round(cy, 8))
            if center_key in seen_centers:
                continue
            seen_centers.add(center_key)

            center_vector = vector_from_tangent_offset(basis, cx, cy)
            candidate = candidate_from_center(center_vector, targets, all_indices, half)
            if anchor not in candidate.covers:
                continue

            candidate_score = score(candidate)
            if candidate_score > best_score:
                best = candidate
                best_score = candidate_score

    return best


def select_priority_pointings(
    targets: list[Target],
    fov_arcmin: float,
    max_pointings: int,
    priority_order: str,
) -> list[Candidate]:
    """Select up to max_pointings by highest-priority uncovered anchor target."""
    reverse = priority_order == "desc"
    order = sorted(
        range(len(targets)),
        key=lambda i: (
            math.isfinite(targets[i].priority),
            targets[i].priority,
        ),
        reverse=reverse,
    )
    order = [i for i in order if math.isfinite(targets[i].priority)]

    selected: list[Candidate] = []
    covered: set[int] = set()

    while len(selected) < max_pointings and len(covered) < len(targets):
        anchor = next((i for i in order if i not in covered), None)
        if anchor is None:
            break

        candidate = best_priority_candidate_for_anchor(
            anchor=anchor,
            targets=targets,
            covered=covered,
            fov_arcmin=fov_arcmin,
        )
        newly_covered = candidate.covers - covered
        selected.append(candidate)
        covered.update(candidate.covers)
        print(
            f"priority pointing {len(selected)}/{max_pointings}: "
            f"anchor={targets[anchor].target_id} "
            f"priority={targets[anchor].priority:.3f} "
            f"newly_covered={len(newly_covered)} "
            f"total_covered={len(covered)}",
            file=sys.stderr,
        )

    return selected


def is_above_priority_threshold(target: Target, threshold: float) -> bool:
    return math.isfinite(target.priority) and target.priority > threshold


def best_threshold_count_candidate_for_anchor(
    anchor: int,
    targets: list[Target],
    covered: set[int],
    fov_arcmin: float,
    priority_threshold: float,
) -> Candidate:
    """
    Choose the best field that includes one above-threshold anchor.

    Score order:
      1. number of not-yet-covered targets with priority > threshold
      2. summed priority of those not-yet-covered above-threshold targets
      3. total number of not-yet-covered targets, including lower-priority extras
      4. edge margin for the anchor
    """
    half = fov_arcmin / 2.0
    if targets[anchor].radius_arcmin > half + 1e-8:
        raise ValueError(
            f"Target {targets[anchor].target_id!r} has radius "
            f"{targets[anchor].radius_arcmin:.6f} arcmin, larger than half the "
            f"{fov_arcmin:.6f} arcmin field. It cannot be fully covered."
        )

    anchor_vector = targets[anchor].vector
    basis = basis_from_vector(anchor_vector)
    local_indices = []
    offsets: dict[int, tuple[float, float]] = {}

    for i, target in enumerate(targets):
        if target.radius_arcmin > half + 1e-8:
            continue
        try:
            x, y = tangent_offsets_arcmin(target.vector, basis)
        except ValueError:
            continue
        max_relevant_offset = math.sqrt(2.0) * fov_arcmin + target.radius_arcmin
        if abs(x) <= max_relevant_offset and abs(y) <= max_relevant_offset:
            local_indices.append(i)
            offsets[i] = (x, y)

    left_edges = sorted(
        set(
            [offsets[i][0] - targets[i].radius_arcmin for i in local_indices]
            + [
                offsets[i][0] + targets[i].radius_arcmin - fov_arcmin
                for i in local_indices
            ]
        )
    )
    bottom_edges = sorted(
        set(
            [offsets[i][1] - targets[i].radius_arcmin for i in local_indices]
            + [
                offsets[i][1] + targets[i].radius_arcmin - fov_arcmin
                for i in local_indices
            ]
        )
    )

    all_indices = list(range(len(targets)))
    best = candidate_from_center(anchor_vector, targets, all_indices, half)

    def score(candidate: Candidate) -> tuple[int, float, int, float]:
        newly_covered = candidate.covers - covered
        high_priority_new = [
            i
            for i in newly_covered
            if is_above_priority_threshold(targets[i], priority_threshold)
        ]
        priority_sum = sum(targets[i].priority for i in high_priority_new)
        return (
            len(high_priority_new),
            priority_sum,
            len(newly_covered),
            min_edge_margin_arcmin(candidate.vector, targets, [anchor], half),
        )

    best_score = score(best)
    seen_centers: set[tuple[float, float]] = set()

    for left in left_edges:
        cx = left + half
        for bottom in bottom_edges:
            cy = bottom + half
            center_key = (round(cx, 8), round(cy, 8))
            if center_key in seen_centers:
                continue
            seen_centers.add(center_key)

            center_vector = vector_from_tangent_offset(basis, cx, cy)
            candidate = candidate_from_center(center_vector, targets, all_indices, half)
            if anchor not in candidate.covers:
                continue

            candidate_score = score(candidate)
            if candidate_score > best_score:
                best = candidate
                best_score = candidate_score

    return best


def select_threshold_count_pointings(
    targets: list[Target],
    fov_arcmin: float,
    max_pointings: int,
    priority_threshold: float,
) -> list[Candidate]:
    """Select up to max_pointings to maximize targets with priority > threshold."""
    eligible = {
        i
        for i, target in enumerate(targets)
        if is_above_priority_threshold(target, priority_threshold)
    }
    if not eligible:
        raise ValueError(
            f"No targets have finite priority greater than {priority_threshold}."
        )

    selected: list[Candidate] = []
    covered: set[int] = set()
    covered_eligible: set[int] = set()

    while len(selected) < max_pointings and covered_eligible != eligible:
        best_anchor = None
        best_candidate = None
        best_score = (-1, -math.inf, -1, -math.inf)

        for anchor in sorted(eligible - covered_eligible, key=lambda i: targets[i].priority, reverse=True):
            candidate = best_threshold_count_candidate_for_anchor(
                anchor=anchor,
                targets=targets,
                covered=covered,
                fov_arcmin=fov_arcmin,
                priority_threshold=priority_threshold,
            )
            newly_covered = candidate.covers - covered
            high_priority_new = [
                i
                for i in newly_covered
                if is_above_priority_threshold(targets[i], priority_threshold)
            ]
            score = (
                len(high_priority_new),
                sum(targets[i].priority for i in high_priority_new),
                len(newly_covered),
                targets[anchor].priority,
            )
            if score > best_score:
                best_anchor = anchor
                best_candidate = candidate
                best_score = score

        if best_candidate is None or best_anchor is None or best_score[0] <= 0:
            break

        newly_covered = best_candidate.covers - covered
        newly_covered_eligible = {
            i
            for i in newly_covered
            if is_above_priority_threshold(targets[i], priority_threshold)
        }
        selected.append(best_candidate)
        covered.update(best_candidate.covers)
        covered_eligible.update(newly_covered_eligible)
        print(
            f"count-threshold pointing {len(selected)}/{max_pointings}: "
            f"anchor={targets[best_anchor].target_id} "
            f"priority={targets[best_anchor].priority:.3f} "
            f"new_high_priority={len(newly_covered_eligible)} "
            f"newly_covered={len(newly_covered)} "
            f"high_priority_covered={len(covered_eligible)}/{len(eligible)}",
            file=sys.stderr,
        )

    return selected


def min_edge_margin_arcmin(
    center_vector: tuple[float, float, float],
    targets: list[Target],
    indices: Iterable[int],
    half_size_arcmin: float,
) -> float:
    basis = basis_from_vector(center_vector)
    margins = []
    for i in indices:
        try:
            x, y = tangent_offsets_arcmin(targets[i].vector, basis)
        except ValueError:
            return -math.inf
        margins.append(
            min(
                half_size_arcmin - abs(x) - targets[i].radius_arcmin,
                half_size_arcmin - abs(y) - targets[i].radius_arcmin,
            )
        )
    return min(margins) if margins else math.inf


def polish_pointing(
    pointing: Candidate,
    targets: list[Target],
    fov_arcmin: float,
) -> Candidate:
    """
    Recenter one selected pointing.

    The pointing's original covered galaxies are treated as required. Among
    centers that still cover those galaxies, prefer centers that cover more
    total catalog sources. If there are only one or two required galaxies, the
    spherical average of their RA/Dec is used unless moving the field can pick
    up additional sources.
    """
    required = set(pointing.covers)
    if not required:
        return pointing

    half = fov_arcmin / 2.0
    all_indices = list(range(len(targets)))
    origin = mean_origin(targets, required)
    basis = basis_from_vector(origin)

    def valid(candidate: Candidate) -> bool:
        return required.issubset(candidate.covers)

    def score(candidate: Candidate) -> tuple[int, int, float]:
        return (
            len(candidate.covers),
            len(candidate.covers - required),
            min_edge_margin_arcmin(candidate.vector, targets, required, half),
        )

    if len(required) <= 2:
        average_vector = origin
        best = candidate_from_center(average_vector, targets, all_indices, half)
        if not valid(best):
            best = pointing
    else:
        best = pointing

    local_indices = []
    offsets: dict[int, tuple[float, float]] = {}
    for i in all_indices:
        if targets[i].radius_arcmin > half + 1e-8:
            continue
        try:
            x, y = tangent_offsets_arcmin(targets[i].vector, basis)
        except ValueError:
            continue
        max_relevant_offset = math.sqrt(2.0) * fov_arcmin + targets[i].radius_arcmin
        if abs(x) <= max_relevant_offset and abs(y) <= max_relevant_offset:
            local_indices.append(i)
            offsets[i] = (x, y)

    for i in required:
        if i not in offsets:
            return best

    left_edges = sorted(
        set(
            [offsets[i][0] - targets[i].radius_arcmin for i in local_indices]
            + [
                offsets[i][0] + targets[i].radius_arcmin - fov_arcmin
                for i in local_indices
            ]
        )
    )
    bottom_edges = sorted(
        set(
            [offsets[i][1] - targets[i].radius_arcmin for i in local_indices]
            + [
                offsets[i][1] + targets[i].radius_arcmin - fov_arcmin
                for i in local_indices
            ]
        )
    )

    best_score = score(best) if valid(best) else (-1, -1, -math.inf)
    prefer_average_center = len(required) <= 2 and len(best.covers) == len(required)
    seen_centers: set[tuple[float, float]] = set()

    for left in left_edges:
        cx = left + half
        for bottom in bottom_edges:
            cy = bottom + half
            center_key = (round(cx, 8), round(cy, 8))
            if center_key in seen_centers:
                continue
            seen_centers.add(center_key)

            center_vector = vector_from_tangent_offset(basis, cx, cy)
            candidate = candidate_from_center(center_vector, targets, all_indices, half)
            if not valid(candidate):
                continue

            candidate_score = score(candidate)
            if candidate_score > best_score:
                if prefer_average_center and len(candidate.covers) <= len(best.covers):
                    continue
                best = candidate
                best_score = candidate_score
                prefer_average_center = False

    return best


def polish_pointings(
    selected: list[Candidate],
    targets: list[Target],
    fov_arcmin: float,
) -> list[Candidate]:
    polished = []
    gained_extra_coverages = 0

    for pointing in selected:
        improved = polish_pointing(pointing, targets, fov_arcmin)
        gained_extra_coverages += max(0, len(improved.covers) - len(pointing.covers))
        polished.append(improved)

    if gained_extra_coverages:
        print(
            f"polish: added {gained_extra_coverages} extra target-field coverages",
            file=sys.stderr,
        )

    return polished


def optimize_pointings(
    targets: list[Target],
    fov_arcmin: float,
    time_limit: float | None,
) -> list[Candidate]:
    selected: list[Candidate] = []
    components = connected_components(targets, fov_arcmin)

    for k, component in enumerate(components, start=1):
        candidates = generate_candidates_for_component(targets, component, fov_arcmin)
        chosen = solve_component(component, candidates, time_limit)
        selected.extend(candidates[j] for j in chosen)
        print(
            f"component {k}/{len(components)}: "
            f"{len(component)} targets, {len(candidates)} candidates, "
            f"{len(chosen)} selected",
            file=sys.stderr,
        )

    return polish_pointings(selected, targets, fov_arcmin)



def write_outputs(
    targets: list[Target],
    selected: list[Candidate],
    fov_arcmin: float,
    pointings_out: str,
    assignments_out: str,
) -> None:
    half = fov_arcmin / 2.0

    with open(pointings_out, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "pointing_id",
                "ra_deg",
                "dec_deg",
                "n_targets_covered",
                "max_priority_covered",
                "sum_priority_covered",
                "covered_target_ids",
            ],
        )
        writer.writeheader()
        for pointing_id, candidate in enumerate(selected, start=1):
            covered_ids = [
                targets[i].target_id for i in sorted(candidate.covers)
            ]
            priorities = [
                targets[i].priority
                for i in candidate.covers
                if math.isfinite(targets[i].priority)
            ]
            writer.writerow(
                {
                    "pointing_id": pointing_id,
                    "ra_deg": f"{candidate.ra_deg:.10f}",
                    "dec_deg": f"{candidate.dec_deg:.10f}",
                    "n_targets_covered": len(candidate.covers),
                    "max_priority_covered": f"{max(priorities):.6f}" if priorities else "",
                    "sum_priority_covered": f"{sum(priorities):.6f}" if priorities else "",
                    "covered_target_ids": ";".join(covered_ids),
                }
            )

    with open(assignments_out, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "target_id",
                "target_ra_deg",
                "target_dec_deg",
                "radius_arcmin",
                "priority",
                "covered_by_selected_pointing",
                "pointing_id",
                "offset_east_arcmin",
                "offset_north_arcmin",
                "edge_margin_arcmin",
            ],
        )
        writer.writeheader()
        for target in targets:
            possible = []
            for pointing_id, candidate in enumerate(selected, start=1):
                if target.index not in candidate.covers:
                    continue
                basis = basis_from_vector(candidate.vector)
                x, y = tangent_offsets_arcmin(target.vector, basis)
                edge_margin = min(
                    half - abs(x) - target.radius_arcmin,
                    half - abs(y) - target.radius_arcmin,
                )
                possible.append((edge_margin, pointing_id, x, y))

            if not possible:
                writer.writerow(
                    {
                        "target_id": target.target_id,
                        "target_ra_deg": f"{target.ra_deg:.10f}",
                        "target_dec_deg": f"{target.dec_deg:.10f}",
                        "radius_arcmin": f"{target.radius_arcmin:.6f}",
                        "priority": f"{target.priority:.6f}",
                        "covered_by_selected_pointing": False,
                        "pointing_id": "",
                        "offset_east_arcmin": "",
                        "offset_north_arcmin": "",
                        "edge_margin_arcmin": "",
                    }
                )
                continue

            edge_margin, pointing_id, x, y = max(possible)
            if edge_margin < -1e-6:
                raise RuntimeError(
                    f"Internal error: target {target.target_id} assignment is outside FOV."
                )

            writer.writerow(
                {
                    "target_id": target.target_id,
                    "target_ra_deg": f"{target.ra_deg:.10f}",
                    "target_dec_deg": f"{target.dec_deg:.10f}",
                    "radius_arcmin": f"{target.radius_arcmin:.6f}",
                    "priority": f"{target.priority:.6f}",
                    "covered_by_selected_pointing": True,
                    "pointing_id": pointing_id,
                    "offset_east_arcmin": f"{x:.6f}",
                    "offset_north_arcmin": f"{y:.6f}",
                    "edge_margin_arcmin": f"{edge_margin:.6f}",
                }
            )


def main() -> list[Candidate]:
    args = parse_args()
    targets = load_catalog(
        args.catalog,
        args.ra_col,
        args.dec_col,
        args.radius_col,
        args.radius_unit,
        args.priority_col,
        args.id_col,
    )
    if args.max_pointings is not None:
        if args.max_pointings < 1:
            raise ValueError("--max-pointings must be at least 1.")
        if args.priority_col is None:
            raise ValueError("--priority-col is required when using --max-pointings.")
        if args.selection_mode == "anchor-priority":
            selected = select_priority_pointings(
                targets=targets,
                fov_arcmin=args.fov_arcmin,
                max_pointings=args.max_pointings,
                priority_order=args.priority_order,
            )
        elif args.selection_mode == "max-count-above-threshold":
            if args.priority_threshold is None:
                raise ValueError(
                    "--priority-threshold is required with "
                    "--selection-mode max-count-above-threshold."
                )
            selected = select_threshold_count_pointings(
                targets=targets,
                fov_arcmin=args.fov_arcmin,
                max_pointings=args.max_pointings,
                priority_threshold=args.priority_threshold,
            )
        else:
            raise ValueError(f"Unsupported selection mode: {args.selection_mode}")
    else:
        selected = optimize_pointings(targets, args.fov_arcmin, args.time_limit)
    write_outputs(
        targets=targets,
        selected=selected,
        fov_arcmin=args.fov_arcmin,
        pointings_out=args.pointings_out,
        assignments_out=args.assignments_out,
    )
    print(f"selected_pointings={len(selected)}")
    print(f"pointings_out={args.pointings_out}")
    print(f"assignments_out={args.assignments_out}")
    return selected

if __name__ == "__main__":
    selected = main()
