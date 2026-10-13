#!/usr/bin/env python

import re
import argparse
from astropy.table import Table
from astropy.coordinates import SkyCoord


def parse_sofia_name(name):
    """
    Parse names like:
        SoFiA J-113830.56+082633.0

    RA part:  hhmmss.ss
    Dec part: +/-ddmmss.s
    """

    name = str(name).strip()


    m = re.search(
        r"J-(\d{1,2})(\d{2})(\d{2}(?:\.\d+)?)([+-])(\d{2})(\d{2})(\d{2}(?:\.\d+)?)",
        name,
    )

    if m is None:
        raise ValueError(f"Could not parse name: {name}")

    rah, ram, ras, decsign, decd, decm, decs = m.groups()

    ra_str = f"{rah}h{ram}m{ras}s"
    dec_str = f"{decsign}{decd}d{decm}m{decs}s"

    c = SkyCoord(ra_str, dec_str, frame="icrs")

    return c.ra.deg, c.dec.deg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument("outfile")
    parser.add_argument("--namecol", default="name")
    parser.add_argument("--informat", default="ascii.basic")
    parser.add_argument("--outformat", default="ascii.commented_header")
    args = parser.parse_args()

    tab = Table.read(
        args.infile,
        format=args.informat,
        guess=False,
        fast_reader=False,
    )

    # just keep galaxies

    # ra_from_name = []
    # dec_from_name = []

    # for name in tab[args.namecol]:
    #     ra, dec = parse_sofia_name(name)
    #     ra_from_name.append(ra)
    #     dec_from_name.append(dec)

    # tab["RA"] = ra_from_name
    # tab["DEC"] = dec_from_name

    tab.write(args.outfile, format=args.outformat, overwrite=True)


if __name__ == "__main__":
    main()
