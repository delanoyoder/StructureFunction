import argparse
import os
import sys

from fits import Fits
from src.plot import Plot
from src.structure_function import RSF, SF
from PIL import Image
import numpy as np


def find_file(file_path):
    if not os.path.isfile(file_path):
        # Try to find the file in the neighboring data directory
        base_dir = os.path.dirname(__file__)
        new_path = os.path.join(base_dir, "..", "data", os.path.basename(file_path))
        if os.path.isfile(new_path):
            print(f"File found in data directory: {new_path}")
            return new_path
        else:
            print(f"File not found at original path or in data directory.")
            sys.exit(1)
    else:
        return file_path


def main(args):
    if args.fits:
        fits_path = find_file(args.fits)
        fits = Fits(fits_path)
        if args.mode == "SF":
            sf = SF(fits.data, args.bins)
            Plot.plot_image(fits)
            Plot.plot_sf(sf)
        elif args.mode == "RSF":
            rsf = RSF(fits.data, args.kernel_size, args.step_size, args.bins)
            Plot.plot_image(fits)
            Plot.plot_rsf(rsf)
    if args.image:
        image_path = find_file(args.image)
        image = Image.open(image_path)
        image = image.convert("RGB")
        image_array = np.array(image)
        if args.channel == "red":
            channel = image_array[:, :, 0]
        elif args.channel == "green":
            channel = image_array[:, :, 1]
        elif args.channel == "blue":
            channel = image_array[:, :, 2]
        if args.mode == "SF":
            sf = SF(channel, args.bins)
            Plot.plot_image(image)
            Plot.plot_sf(sf)
        elif args.mode == "RSF":
            breakpoint()
            rsf = RSF(channel, args.kernel_size, args.step_size, args.bins)
            Plot.plot_image(image)
            Plot.plot_rsf(rsf)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # --fits flag, which requires an argument
    parser.add_argument(
        "--fits", type=str, required=False, help="Path to the fits file: sample.fits"
    )

    # --image flag, which requires an argument
    parser.add_argument(
        "--image", type=str, required=False, help="Path to the image file: sample.png"
    )

    # --bins flag, which requires an argument
    parser.add_argument(
        "--bins", type=int, required=False, help="Number of bins to use: 10, 20, 30"
    )

    # --mode flag, which requires an argument
    parser.add_argument(
        "--mode",
        type=str,
        required=False,
        help="Type of structure function to use: SF, RSF, SCASF, SCARSF",
    )

    # --channel flag, which requires an argument
    parser.add_argument(
        "--channel",
        type=str,
        default="red",
        required=False,
        help="Channel to use: red, green, blue",
    )

    parser.add_argument(
        "--kernel-size",
        type=int,
        default=5,
        required=False,
        help="Size of kernel: 5, 10, 15",
    )

    parser.add_argument(
        "--step-size", type=int, default=1, required=False, help="Size of step: 1, 2, 3"
    )

    args = parser.parse_args()

    # If neither fits nor image path is provided, raise an error and exit
    if args.fits is None and args.image is None:
        print("Error: At least one of --fits or --image is required.")
        sys.exit(1)

    main(args)
