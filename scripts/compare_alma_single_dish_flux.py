#!/usr/bin/env python
"""Compare ALMA CO(1-0) moment-0 central fluxes to single-dish beams.

The ALMA maps are expected to have moment-0 units of Jy km/s / beam and
beam keywords BMAJ/BMIN in degrees.  The script smooths each image to a
requested single-dish HPBW, then reports the value at the central pixel in
Jy km/s per single-dish beam.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import gaussian_filter


DEFAULT_ALMA_DIR = Path("/Users/rfinn/research/Virgo/alma/combes_data")
DEFAULT_RESULTS_DIR = DEFAULT_ALMA_DIR.parent / "alma_single_dish_comparison"
DEFAULT_OUTPUT = DEFAULT_RESULTS_DIR / "alma_center_fluxes.csv"
DEFAULT_COMPARISON = Path(__file__).with_name("alma_single_dish_measurements.csv")
FWHM_TO_SIGMA = 1.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))


@dataclass(frozen=True)
class Target:
    alma_image: str
    single_dish_observation: str
    beam_key: str


DEFAULT_TARGETS = (
    Target("n5348-co10-mean.fits", "IRAM_30m", "iram30m"),
    Target("n5356-co10-mean.fits", "KPNO_12m", "kpno12m"),
    Target("n5470-co10-mean.fits", "IRAM_30m", "iram30m"),
    Target("n5560-co10-mean.fits", "KPNO_12m", "kpno12m"),
    Target("n5566-co10-mean.fits", "IRAM_30m", "iram30m"),
    Target("n5577-co10-mean.fits", "KPNO_12m", "kpno12m"),
    Target("u9661-co10-mean.fits", "IRAM_30m", "iram30m"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Smooth ALMA moment-zero maps to a single-dish HPBW and measure "
            "the central-pixel flux."
        )
    )
    parser.add_argument(
        "--alma-dir",
        type=Path,
        default=DEFAULT_ALMA_DIR,
        help=f"Directory containing ALMA FITS files. Default: {DEFAULT_ALMA_DIR}",
    )
    parser.add_argument(
        "--images",
        nargs="*",
        help=(
            "Optional image filenames to process. Defaults to the seven ALMA/"
            "single-dish pairs listed in the script."
        ),
    )
    parser.add_argument(
        "--single-dish-beam-arcsec",
        type=float,
        help=(
            "Use this single-dish HPBW for all images, in arcsec. If omitted, "
            "the built-in IRAM/KPNO defaults are used for the default image list."
        ),
    )
    parser.add_argument(
        "--iram30m-beam-arcsec",
        type=float,
        default=21.0,
        help="IRAM 30m CO(1-0) HPBW in arcsec. Default: 21",
    )
    parser.add_argument(
        "--kpno12m-beam-arcsec",
        type=float,
        default=55.0,
        help="KPNO 12m CO(1-0) HPBW in arcsec. Default: 55",
    )
    parser.add_argument(
        "--center",
        choices=("image", "crpix"),
        default="image",
        help="Central pixel definition. Default: image center.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"CSV output path. Default: {DEFAULT_OUTPUT}",
    )
    parser.add_argument(
        "--comparison-table",
        type=Path,
        default=DEFAULT_COMPARISON,
        help=(
            "CSV/ECSV/FITS table of previous single-dish fluxes. Defaults to "
            "the supplied measurements in alma_single_dish_measurements.csv. "
            "Rows are matched by alma_image, image, filename, source_name, "
            "galaxy, or name when available."
        ),
    )
    parser.add_argument(
        "--comparison-error-column",
        default="errSCO_10",
        help="Optional flux uncertainty column. Default: errSCO_10",
    )
    parser.add_argument(
        "--comparison-flux-column",
        default="SCO_10",
        help="Flux column to read from --comparison-table. Default: SCO_10",
    )
    parser.add_argument(
        "--write-convolved",
        action="store_true",
        help="Write smoothed, target-beam-scaled FITS images.",
    )
    parser.add_argument(
        "--convolved-dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR / "convolved",
        help=f"Directory for --write-convolved outputs. Default: {DEFAULT_RESULTS_DIR / 'convolved'}",
    )
    parser.add_argument(
        "--gaussian-truncate",
        type=float,
        default=4.0,
        help="Gaussian kernel radius in sigma. Default: 4",
    )
    return parser.parse_args()


def target_rows(args: argparse.Namespace) -> list[dict[str, object]]:
    beam_by_key = {
        "iram30m": args.iram30m_beam_arcsec,
        "kpno12m": args.kpno12m_beam_arcsec,
    }
    target_by_image = {target.alma_image: target for target in DEFAULT_TARGETS}

    if args.images:
        rows = []
        for image in args.images:
            target = target_by_image.get(Path(image).name)
            if target is None:
                if args.single_dish_beam_arcsec is None:
                    raise ValueError(
                        "--single-dish-beam-arcsec is required for images "
                        "outside the built-in target list."
                    )
                rows.append(
                    {
                        "alma_image": image,
                        "single_dish_observation": "custom",
                        "single_dish_beam_arcsec": args.single_dish_beam_arcsec,
                    }
                )
            else:
                rows.append(
                    {
                        "alma_image": target.alma_image,
                        "single_dish_observation": target.single_dish_observation,
                        "single_dish_beam_arcsec": args.single_dish_beam_arcsec
                        or beam_by_key[target.beam_key],
                    }
                )
        return rows

    return [
        {
            "alma_image": target.alma_image,
            "single_dish_observation": target.single_dish_observation,
            "single_dish_beam_arcsec": args.single_dish_beam_arcsec
            or beam_by_key[target.beam_key],
        }
        for target in DEFAULT_TARGETS
    ]


def read_image_plane(path: Path) -> tuple[np.ndarray, fits.Header, tuple[int, ...]]:
    with fits.open(path) as hdul:
        data = np.asarray(hdul[0].data, dtype=float)
        header = hdul[0].header.copy()

    original_shape = data.shape
    plane = np.squeeze(data)
    if plane.ndim != 2:
        raise ValueError(f"{path} does not reduce to one 2D image plane: {original_shape}")
    return plane, header, original_shape


def celestial_pixel_scale_arcsec(header: fits.Header) -> float:
    wcs = WCS(header).celestial
    scales_deg = np.asarray(proj_plane_pixel_scales(wcs), dtype=float)
    if len(scales_deg) < 2 or not np.all(np.isfinite(scales_deg[:2])):
        cdelt1 = abs(float(header["CDELT1"]))
        cdelt2 = abs(float(header["CDELT2"]))
        scales_deg = np.array([cdelt1, cdelt2])
    return float(np.sqrt(np.prod(scales_deg[:2])) * 3600.0)


def center_pixel(header: fits.Header, shape: tuple[int, int], mode: str) -> tuple[int, int]:
    ny, nx = shape
    if mode == "crpix" and "CRPIX1" in header and "CRPIX2" in header:
        x = int(round(float(header["CRPIX1"]))) - 1
        y = int(round(float(header["CRPIX2"]))) - 1
    else:
        x = nx // 2
        y = ny // 2
    return min(max(x, 0), nx - 1), min(max(y, 0), ny - 1)


def nan_gaussian_filter(data: np.ndarray, sigma_pix: float, truncate: float) -> np.ndarray:
    finite = np.isfinite(data)
    filled = np.where(finite, data, 0.0)
    weights = finite.astype(float)

    smoothed = gaussian_filter(
        filled,
        sigma=sigma_pix,
        mode="constant",
        cval=0.0,
        truncate=truncate,
    )
    smoothed_weights = gaussian_filter(
        weights,
        sigma=sigma_pix,
        mode="constant",
        cval=0.0,
        truncate=truncate,
    )

    output = np.full_like(smoothed, np.nan, dtype=float)
    np.divide(
        smoothed,
        smoothed_weights,
        out=output,
        where=smoothed_weights > 0,
    )
    return output


def restored_shape(plane: np.ndarray, original_shape: tuple[int, ...]) -> np.ndarray:
    if len(original_shape) == 2:
        return plane
    if len(original_shape) == 3 and original_shape[0] == 1:
        return plane[np.newaxis, :, :]
    return plane.reshape(original_shape)


def write_convolved_fits(
    output_path: Path,
    smoothed_target_beam: np.ndarray,
    original_shape: tuple[int, ...],
    header: fits.Header,
    target_beam_arcsec: float,
    kernel_fwhm_arcsec: float,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_header = header.copy()
    out_header["BMAJ"] = (target_beam_arcsec / 3600.0, "Target circular beam FWHM [deg]")
    out_header["BMIN"] = (target_beam_arcsec / 3600.0, "Target circular beam FWHM [deg]")
    out_header["BPA"] = (0.0, "Target circular beam position angle [deg]")
    out_header["HISTORY"] = (
        f"Smoothed to {target_beam_arcsec:.6g} arcsec HPBW with "
        f"{kernel_fwhm_arcsec:.6g} arcsec FWHM Gaussian kernel."
    )
    out_header["HISTORY"] = (
        "Scaled by target/original beam area; units are per target beam."
    )
    fits.writeto(
        output_path,
        restored_shape(smoothed_target_beam, original_shape),
        out_header,
        overwrite=True,
    )


def normalize_key(value: object) -> str:
    key = Path(str(value)).name.lower()
    for token in ("-co10-mean.fits", ".fits", "ngc", "ugc", " "):
        key = key.replace(token, "")
    if len(key) > 1 and key[0] in ("n", "u") and key[1:].isdigit():
        key = key[1:]
    return key.strip()


def load_comparison_fluxes(
    path: Path, flux_column: str, error_column: str
) -> dict[str, dict[str, object]]:
    table = Table.read(path)
    if flux_column not in table.colnames:
        raise ValueError(
            f"Column {flux_column!r} not found in {path}. "
            f"Available columns: {', '.join(table.colnames)}"
        )

    key_columns = [
        "alma_image",
        "image",
        "filename",
        "source_name",
        "galaxy",
        "name",
        "NEDname",
    ]
    available_key_columns = [col for col in key_columns if col in table.colnames]
    if not available_key_columns:
        raise ValueError(
            f"No recognized key columns found in {path}. Expected one of: "
            f"{', '.join(key_columns)}"
        )

    comparison = {}
    for row in table:
        flux = row[flux_column]
        try:
            flux = float(flux)
        except (TypeError, ValueError):
            continue
        if not np.isfinite(flux):
            continue
        error = np.nan
        if error_column in table.colnames:
            try:
                error = float(row[error_column])
            except (TypeError, ValueError):
                pass
        measurement = {
            "flux": flux,
            "error": error,
            "pointing_mode": str(row["pointing_mode"])
            if "pointing_mode" in table.colnames else "unknown",
            "telescope": str(row["single_dish_observation"])
            if "single_dish_observation" in table.colnames else "unknown",
            "beam_arcsec": float(row["single_dish_beam_arcsec"])
            if "single_dish_beam_arcsec" in table.colnames else np.nan,
        }
        for col in available_key_columns:
            comparison[normalize_key(row[col])] = measurement
    return comparison


def measure_target(
    row: dict[str, object],
    args: argparse.Namespace,
    comparison_fluxes: dict[str, dict[str, object]] | None,
) -> dict[str, object]:
    image_name = str(row["alma_image"])
    path = Path(image_name)
    if not path.is_absolute():
        path = args.alma_dir / path

    data, header, original_shape = read_image_plane(path)
    bmaj_arcsec = float(header["BMAJ"]) * 3600.0
    bmin_arcsec = float(header["BMIN"]) * 3600.0
    alma_beam_arcsec = math.sqrt(bmaj_arcsec * bmin_arcsec)
    target_beam_arcsec = float(row["single_dish_beam_arcsec"])
    if target_beam_arcsec <= alma_beam_arcsec:
        raise ValueError(
            f"{path.name}: target beam ({target_beam_arcsec:.3f} arcsec) must be "
            f"larger than ALMA geometric beam ({alma_beam_arcsec:.3f} arcsec)."
        )

    kernel_fwhm_arcsec = math.sqrt(target_beam_arcsec**2 - alma_beam_arcsec**2)
    pixel_scale_arcsec = celestial_pixel_scale_arcsec(header)
    kernel_sigma_pix = kernel_fwhm_arcsec * FWHM_TO_SIGMA / pixel_scale_arcsec

    smoothed_alma_beam = nan_gaussian_filter(
        data,
        sigma_pix=kernel_sigma_pix,
        truncate=args.gaussian_truncate,
    )
    beam_area_ratio = target_beam_arcsec**2 / (bmaj_arcsec * bmin_arcsec)
    smoothed_target_beam = smoothed_alma_beam * beam_area_ratio

    x_center, y_center = center_pixel(header, data.shape, args.center)
    raw_center = float(data[y_center, x_center])
    smooth_center_alma_beam = float(smoothed_alma_beam[y_center, x_center])
    smooth_center_target_beam = float(smoothed_target_beam[y_center, x_center])

    previous_flux = np.nan
    previous_error = np.nan
    pointing_mode = "unknown"
    comparison_note = "no matching measurement"
    delta_flux = np.nan
    flux_ratio = np.nan
    if comparison_fluxes is not None:
        for key in (path.name, path.stem, image_name):
            measurement = comparison_fluxes.get(normalize_key(key))
            if measurement is not None:
                previous_flux = measurement["flux"]
                previous_error = measurement["error"]
                pointing_mode = measurement["pointing_mode"]
                comparison_note = {
                    "single": "central beam vs single pointing",
                    "multiple": "central beam vs multiple pointings; not directly comparable",
                }.get(pointing_mode, "pointing geometry unspecified")
                measured_beam = measurement["beam_arcsec"]
                if np.isfinite(measured_beam) and not np.isclose(
                    target_beam_arcsec, measured_beam
                ):
                    comparison_note += "; prediction beam differs from measured beam"
                if measurement["telescope"] not in (
                    "unknown", row["single_dish_observation"]
                ):
                    comparison_note += "; measurement telescope differs from prediction"
                break
        if np.isfinite(previous_flux):
            delta_flux = smooth_center_target_beam - previous_flux
            if previous_flux != 0:
                flux_ratio = smooth_center_target_beam / previous_flux

    if args.write_convolved:
        out_name = f"{path.stem}_{target_beam_arcsec:g}arcsec.fits"
        write_convolved_fits(
            args.convolved_dir / out_name,
            smoothed_target_beam,
            original_shape,
            header,
            target_beam_arcsec,
            kernel_fwhm_arcsec,
        )

    return {
        "alma_image": path.name,
        "single_dish_observation": row["single_dish_observation"],
        "single_dish_beam_arcsec": target_beam_arcsec,
        "header_BMAJ_deg": float(header["BMAJ"]),
        "header_BMIN_deg": float(header["BMIN"]),
        "alma_bmaj_arcsec": bmaj_arcsec,
        "alma_bmin_arcsec": bmin_arcsec,
        "alma_geometric_beam_arcsec": alma_beam_arcsec,
        "kernel_fwhm_arcsec": kernel_fwhm_arcsec,
        "pixel_scale_arcsec": pixel_scale_arcsec,
        "kernel_sigma_pix": kernel_sigma_pix,
        "center_x_0based": x_center,
        "center_y_0based": y_center,
        "center_x_fits": x_center + 1,
        "center_y_fits": y_center + 1,
        "raw_center_jy_km_s_per_alma_beam": raw_center,
        "smoothed_center_jy_km_s_per_alma_beam": smooth_center_alma_beam,
        "beam_area_ratio_target_over_alma": beam_area_ratio,
        "alma_center_jy_km_s_per_single_dish_beam": smooth_center_target_beam,
        "previous_single_dish_flux": previous_flux,
        "previous_single_dish_flux_error": previous_error,
        "single_dish_pointing_mode": pointing_mode,
        "comparison_note": comparison_note,
        "alma_minus_previous_flux": delta_flux,
        "alma_over_previous_flux": flux_ratio,
        "bunit": header.get("BUNIT", ""),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    rows = target_rows(args)
    comparison_fluxes = None
    if args.comparison_table:
        comparison_fluxes = load_comparison_fluxes(
            args.comparison_table,
            args.comparison_flux_column,
            args.comparison_error_column,
        )

    measurements = [
        measure_target(row, args, comparison_fluxes=comparison_fluxes)
        for row in rows
    ]
    write_csv(args.output, measurements)

    display_columns = [
        "alma_image",
        "single_dish_observation",
        "single_dish_beam_arcsec",
        "alma_geometric_beam_arcsec",
        "kernel_fwhm_arcsec",
        "alma_center_jy_km_s_per_single_dish_beam",
        "previous_single_dish_flux",
        "previous_single_dish_flux_error",
        "single_dish_pointing_mode",
        "alma_over_previous_flux",
    ]
    Table(rows=measurements)[display_columns].pprint(max_width=-1, max_lines=-1)
    print("\nFluxes and measurement uncertainties: Jy km/s. Predictions: per target beam.")
    print("KPNO multiple-pointing ratios are not direct single-pointing comparisons.")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
