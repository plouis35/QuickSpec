from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
from astropy.nddata import CCDData

from img.img_utils import ImagesCombiner
from tools.crop_fits import _output_path, main


@pytest.mark.parametrize(
    ("source_name", "expected_name"),
    [
        ("image-12.fits", "image-cropped-12.fits"),
        ("target.v2-007.fit", "target.v2-cropped-007.fit"),
        ("image.fits", "image_cropped.fits"),
        ("image-12.fits.gz", "image-cropped-12.fits.gz"),
    ],
)
def test_output_path_inserts_cropped_before_trailing_dash_number(
    source_name, expected_name
):
    assert _output_path(Path(source_name)).name == expected_name


def test_y_crop_keeps_full_width(synthetic_image):
    images = ImagesCombiner([synthetic_image], ["synthetic.fits"])

    result = images.y_crop("0.5,0.4")

    assert result is images
    assert images[0].shape == (40, 200)
    np.testing.assert_array_equal(
        images[0].data,
        synthetic_image.data[30:70, :],
    )


def test_cli_expands_glob_and_writes_cropped_file(
    image_fits_file, synthetic_image, capsys
):
    source = Path(image_fits_file)

    arguments = [
        str(source.parent / "test_*.fits"),
        "--center", "0.5",
        "--size", "0.2",
    ]
    result = main(arguments)

    output = source.with_name("test_image_cropped.fits")
    assert result == 0
    assert output.exists()
    cropped = CCDData.read(output, unit=u.Unit("adu"))
    assert cropped.shape == (40, 200)
    np.testing.assert_array_equal(
        cropped.data,
        synthetic_image.data[30:70, :],
    )

    with pytest.raises(SystemExit) as error:
        main(arguments)

    assert error.value.code == 2
    assert "output file already exists" in capsys.readouterr().err
