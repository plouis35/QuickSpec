#!/usr/bin/env python3
"""Convert Bayer RAW images to monochrome FITS with 2x2 CFA binning.

Requires rawpy, numpy, and astropy:
    python -m pip install rawpy numpy astropy
"""

import argparse
from pathlib import Path

import numpy as np
import rawpy
from astropy.io import fits


def convert(raw_path, output_path, pixel_size_um=None, overwrite=False):
    with rawpy.imread(str(raw_path)) as raw:
        image = raw.raw_image_visible
        height, width = image.shape
        binned_height = height - height % 2
        binned_width = width - width % 2
        if binned_height == 0 or binned_width == 0:
            raise ValueError(f"RAW image is too small to bin: {raw_path}")

        pattern = raw.raw_pattern
        if pattern is None or pattern.shape != (2, 2):
            raise ValueError(f"Image does not have a supported 2x2 Bayer pattern: {raw_path}")

        black_levels = raw.black_level_per_channel
        top = raw.sizes.top_margin
        left = raw.sizes.left_margin
        binned = np.zeros(
            (binned_height // 2, binned_width // 2), dtype=np.int32
        )
        for dy in range(2):
            for dx in range(2):
                channel = int(pattern[(top + dy) % 2, (left + dx) % 2])
                black = black_levels[channel] if channel < len(black_levels) else 0
                plane = image[dy:binned_height:2, dx:binned_width:2].astype(
                    np.int32
                )
                binned += plane - black

        header = fits.Header()
        header["BUNIT"] = ("ADU", "Signal above RAW black level")
        header["BAYERBIN"] = (2, "2x2 Bayer pixels summed per output pixel")
        header["XBINNING"] = (2, "Horizontal binning factor")
        header["YBINNING"] = (2, "Vertical binning factor")
        header["BLACKSUB"] = (True, "Per-channel RAW black level subtracted")
        header["RAWNX"] = (binned_width, "Visible RAW width used before binning")
        header["RAWNY"] = (binned_height, "Visible RAW height used before binning")
        if binned_width != width:
            header["TRIMX"] = (width - binned_width, "Pixels trimmed from right edge")
        if binned_height != height:
            header["TRIMY"] = (height - binned_height, "Pixels trimmed from bottom edge")

        if pixel_size_um is not None:
            binned_pixel_size = 2 * pixel_size_um
            header["XPIXSZ"] = (binned_pixel_size, "Pixel width in microns")
            header["YPIXSZ"] = (binned_pixel_size, "Pixel height in microns")

        fits.writeto(output_path, binned, header=header, overwrite=overwrite)


def main():
    parser = argparse.ArgumentParser(
        description="Convert RAW Bayer images to monochrome FITS using 2x2 binning."
    )
    parser.add_argument("raw_files", nargs="+", type=Path, help="Input RAW file(s)")
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help="Output directory (defaults to each input file's directory)",
    )
    parser.add_argument(
        "--pixel-size-um",
        type=float,
        help="Physical size of one RAW photosite in microns; FITS pixel size is doubled",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Replace existing FITS files"
    )
    args = parser.parse_args()

    if args.pixel_size_um is not None and args.pixel_size_um <= 0:
        parser.error("--pixel-size-um must be greater than zero")
    if args.output_dir is not None:
        args.output_dir.mkdir(parents=True, exist_ok=True)

    for raw_path in args.raw_files:
        output_dir = args.output_dir or raw_path.parent
        output_path = output_dir / f"{raw_path.stem}.fits"
        if output_path.exists() and not args.overwrite:
            parser.error(f"{output_path} already exists; use --overwrite to replace it")
        convert(raw_path, output_path, args.pixel_size_um, args.overwrite)
        print(f"{raw_path} -> {output_path}")


if __name__ == "__main__":
    main()
