"""
CLI to crop FITS images vertically while keeping their full width.
Center and size are fractions of image height.
"""

import argparse
import glob
import re
from pathlib import Path

from img.img_utils import Images

def _expand_inputs(patterns: list[str]) -> list[Path]:
    """Expand shell-style patterns and keep each input only once."""
    files = dict.fromkeys(
        Path(match)
        for pattern in patterns
        for match in glob.glob(pattern, recursive=True)
    )
    return sorted(files)


def _output_path(source: Path) -> Path:
    """ compute output filenames depending according to sequence nb if any"""
    fits_extensions = (".fits.gz", ".fit.gz", ".fts.gz", ".fits", ".fit", ".fts")
    extension = next(
        (
            source.name[-len(candidate):]
            for candidate in fits_extensions
            if source.name.lower().endswith(candidate)
        ),
        source.suffix,
    )
    stem = source.name[:-len(extension)] if extension else source.name
    match = re.search(r"(-\d+)$", stem)
    if match:
        stem = f"{stem[:match.start()]}-cropped{match.group(1)}"
    else:
        stem = f"{stem}_cropped"
    return source.with_name(f"{stem}{extension}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Crop FITS images vertically while keeping their full width. "
            "Center and size are fractions of image height."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="FITS file paths or wildcard patterns",
    )
    parser.add_argument(
        "--center",
        type=float,
        metavar="Y",
        default=0.5,
        help="relative vertical crop center, from 0 to 1 (default: 0.5)",
    )
    parser.add_argument(
        "--size",
        type=float,
        metavar="RATIO",
        default=0.2,
        help=(
            "fraction of image height kept on either side of the center "
            "(default: 0.2)"
        ),
    )
    args = parser.parse_args(argv)

    if not 0 <= args.center <= 1:
        parser.error("center must be between 0 and 1")
    if not 0 < args.size <= 1:
        parser.error("size must be greater than 0 and at most 1")
    if args.center - args.size < 0 or args.center + args.size > 1:
        parser.error("crop region extends beyond the image; adjust center or size")

    files = _expand_inputs(args.inputs)
    if not files:
        parser.error("no input files matched")

    for source in files:
        output = _output_path(source)
        if output.exists():
            parser.error(f"output file already exists: {output}")
        images = Images.from_fits_by_names_list([str(source)])
        try:
            images.y_crop(f"{args.center},{2 * args.size}")
        except ValueError as error:
            parser.error(f"{source}: {error}")
        images[0].write(output)
        print(f"{source} -> {output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
