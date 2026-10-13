#!/usr/bin/env python
"""
Find a minimum set of square telescope pointings that covers a target catalog.

Example:
    python lbt_square_cover.py vf_v2_main_meerkat.fits --ra RA --dec DEC --id VFID \
        --name objname --fov 25 --outroot lbt_square25

This treats the field as an axis-aligned square in a local tangent-plane
approximation, with side length in arcmin.  It covers target centers only.
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.optimize import milp, LinearConstraint, Bounds


def minimum_square_cover(tab, ra_col="RA", dec_col="DEC", id_col=None, name_col=None, fov_arcmin=25.0):
    ra = np.asarray(tab[ra_col], dtype=float)
    dec = np.asarray(tab[dec_col], dtype=float)
    n = len(ra)

    if id_col is None:
        ids = np.array([str(i) for i in range(n)])
    else:
        ids = np.array([str(x).strip() for x in tab[id_col]])

    if name_col is None:
        names = ids.copy()
    else:
        names = np.array([str(x).strip() for x in tab[name_col]])

    ra0 = float(np.mean(ra))
    dec0 = float(np.mean(dec))

    # Local tangent-plane coordinates in arcmin.
    x = (ra - ra0) * np.cos(np.deg2rad(dec0)) * 60.0
    y = (dec - dec0) * 60.0
    L = float(fov_arcmin)
    eps = 1e-9

    # For axis-aligned fixed-size squares, the set of possible coverages only
    # changes when a square edge passes through a target.  It is therefore
    # sufficient to test left edges at x_i and x_i - L, and bottom edges at
    # y_i and y_i - L.
    lefts = np.unique(np.concatenate([x, x - L]))
    bottoms = np.unique(np.concatenate([y, y - L]))

    unique = {}
    for lx in lefts:
        inx = np.where((x >= lx - eps) & (x <= lx + L + eps))[0]
        if len(inx) == 0:
            continue
        for by in bottoms:
            inds = inx[(y[inx] >= by - eps) & (y[inx] <= by + L + eps)]
            if len(inds) == 0:
                continue
            mask = 0
            for i in inds:
                mask |= (1 << int(i))
            if mask not in unique:
                unique[mask] = (lx + L / 2, by + L / 2, len(inds))

    masks = np.array(list(unique.keys()), dtype=object)
    counts = np.array([v[2] for v in unique.values()])
    m = len(masks)

    A = np.zeros((n, m), dtype=float)
    for j, mask in enumerate(masks):
        for i in range(n):
            if (mask >> i) & 1:
                A[i, j] = 1.0

    result = milp(
        c=np.ones(m),
        constraints=LinearConstraint(A, lb=np.ones(n), ub=np.full(n, np.inf)),
        bounds=Bounds(0, 1),
        integrality=np.ones(m),
        options={"mip_rel_gap": 0},
    )

    if not result.success:
        raise RuntimeError(f"MILP failed: {result.message}")

    selected = np.where(result.x > 0.5)[0]

    # Recenter selected boxes on their covered targets to improve edge margins.
    pointings = []
    for j in selected:
        inds = np.array([i for i in range(n) if (masks[j] >> i) & 1], dtype=int)
        cx = 0.5 * (x[inds].min() + x[inds].max())
        cy = 0.5 * (y[inds].min() + y[inds].max())

        dx = np.abs(x - cx)
        dy = np.abs(y - cy)
        infield = np.where((dx <= L / 2 + eps) & (dy <= L / 2 + eps))[0]
        margin = float(np.min(np.minimum(L / 2 - dx[infield], L / 2 - dy[infield])))

        ra_c = ra0 + cx / (60.0 * np.cos(np.deg2rad(dec0)))
        dec_c = dec0 + cy / 60.0

        pointings.append({
            "x": float(cx),
            "y": float(cy),
            "ra": float(ra_c),
            "dec": float(dec_c),
            "infield": infield,
            "margin": margin,
        })

    # Sort by RA.
    order = np.argsort([p["ra"] for p in pointings])
    pointings = [pointings[i] for i in order]

    # Assign each target to the selected pointing where it has the largest edge margin.
    assignment = []
    for i in range(n):
        best = None
        best_margin = -np.inf
        for k, p in enumerate(pointings):
            dx_i = abs(x[i] - p["x"])
            dy_i = abs(y[i] - p["y"])
            if dx_i <= L / 2 + eps and dy_i <= L / 2 + eps:
                margin_i = min(L / 2 - dx_i, L / 2 - dy_i)
                if margin_i > best_margin:
                    best = k
                    best_margin = margin_i
        if best is None:
            raise RuntimeError(f"Target {ids[i]} is uncovered.")
        assignment.append((best, float(best_margin)))

    pointing_rows = []
    assignment_rows = []
    for k, p in enumerate(pointings, start=1):
        infield = list(map(int, p["infield"]))
        assigned = [i for i, a in enumerate(assignment) if a[0] == k - 1]

        for i in assigned:
            assignment_rows.append({
                "POINTING": k,
                "ID": ids[i],
                "NAME": names[i],
                "RA": ra[i],
                "DEC": dec[i],
                "EDGE_MARGIN_ARCMIN": assignment[i][1],
            })

        pointing_rows.append({
            "POINTING": k,
            "RA_CENTER": p["ra"],
            "DEC_CENTER": p["dec"],
            "N_ASSIGNED": len(assigned),
            "N_IN_FIELD": len(infield),
            "MIN_MARGIN_ASSIGNED_ARCMIN": min([assignment[i][1] for i in assigned]) if assigned else np.nan,
            "MIN_MARGIN_ANY_IN_FIELD_ARCMIN": p["margin"],
            "ASSIGNED_IDS": ",".join(ids[assigned]),
            "ASSIGNED_NAMES": ",".join(names[assigned]),
            "ALL_IN_FIELD_IDS": ",".join(ids[infield]),
            "ALL_IN_FIELD_NAMES": ",".join(names[infield]),
        })

    return pd.DataFrame(pointing_rows), pd.DataFrame(assignment_rows), result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("catalog")
    parser.add_argument("--ra", default="RA")
    parser.add_argument("--dec", default="DEC")
    parser.add_argument("--id", default=None)
    parser.add_argument("--name", default=None)
    parser.add_argument("--fov", type=float, default=25.0)
    parser.add_argument("--outroot", default="square_cover")
    args = parser.parse_args()

    tab = Table.read(args.catalog)
    pointings, assignments, result = minimum_square_cover(
        tab,
        ra_col=args.ra,
        dec_col=args.dec,
        id_col=args.id,
        name_col=args.name,
        fov_arcmin=args.fov,
    )

    pointings.to_csv(f"{args.outroot}_pointings.csv", index=False)
    assignments.to_csv(f"{args.outroot}_assignments.csv", index=False)

    print(f"Optimal number of pointings: {len(pointings)}")
    print(f"Wrote {args.outroot}_pointings.csv")
    print(f"Wrote {args.outroot}_assignments.csv")
    return pointings, assignments


if __name__ == "__main__":
    #pointints, assignments = main()
    parser = argparse.ArgumentParser()
    parser.add_argument("catalog")
    parser.add_argument("--ra", default="RA")
    parser.add_argument("--dec", default="DEC")
    parser.add_argument("--id", default=None)
    parser.add_argument("--name", default=None)
    parser.add_argument("--fov", type=float, default=25.0)
    parser.add_argument("--outroot", default="square_cover")
    args = parser.parse_args()

    tab = Table.read(args.catalog)
    pointings, assignments, result = minimum_square_cover(
        tab,
        ra_col=args.ra,
        dec_col=args.dec,
        id_col=args.id,
        name_col=args.name,
        fov_arcmin=args.fov,
    )

    pointings.to_csv(f"{args.outroot}_pointings.csv", index=False)
    assignments.to_csv(f"{args.outroot}_assignments.csv", index=False)

    print(f"Optimal number of pointings: {len(pointings)}")
    print(f"Wrote {args.outroot}_pointings.csv")
    print(f"Wrote {args.outroot}_assignments.csv")
