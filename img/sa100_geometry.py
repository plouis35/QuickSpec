"""
sa100_geometry — StarAnalyzer slitless spectrograph processing

Given two user-clicked points (order 0 and a point on the spectrum),
this module computes:
  - the dispersion angle
  - the rectified (rotated) image aligned with the dispersion axis
  - the spectrum and sky crop zones
  - the linear wavelength calibration axis

SA100 dispersion formula (linear approximation):
    dispersion [Å/px] = (pixel_size_µm / distance_mm) × K
where K ≈ 5765 is derived from the grating equation at 5500 Å.
This is accurate to ~2% across the visible range, sufficient for
the SA100's ~100 Å/px resolution.
"""
import logging
import math
import numpy as np
from datetime import datetime

from astropy import units as u
from astropy.nddata import CCDData
from astropy.table import Table
from specutils import Spectrum
from app.config import Config
from skimage.transform import rotate as ski_rotate

# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def compute_angle(x0: float, y0: float, x1: float, y1: float) -> float:
    """
    Compute the dispersion angle in degrees from two clicked points.

    Args:
        x0, y0: pixel coordinates of order 0
        x1, y1: pixel coordinates of a point on the spectrum (order 1)

    Returns:
        float: angle in degrees (positive = counterclockwise)
    """
    dx = x1 - x0
    dy = y1 - y0
    angle = math.degrees(math.atan2(dy, dx))
    logging.info(f"SA100: dispersion angle = {angle:.2f}°  (dx={dx:.1f}, dy={dy:.1f})")
    return angle


def rectify(image: CCDData, angle: float) -> CCDData:
    """
    Rotate the image so the dispersion axis is horizontal.

    Uses skimage.transform.rotate with bicubic interpolation and
    constant (zero) border fill — no cyclic artefacts.

    Args:
        image (CCDData): raw 2D image
        angle (float):  dispersion angle in degrees (from compute_angle)

    Returns:
        CCDData: rotated image, same shape as input
    """
    rotated = ski_rotate(
        image.data.astype(np.float32),
#        angle=-angle,          # counter-rotate to align spectrum horizontally
        angle=angle,          # counter-rotate to align spectrum horizontally
        resize=False,          # keep original shape
        order=3,               # bicubic
        mode='constant',
        cval=0.0,
        preserve_range=True,
    )
    result = CCDData(rotated, unit=image.unit, header=image.header.copy())
    result.header['SA100ROT'] = round(angle, 3)
    logging.info(f"SA100: image rectified (rotation = {-angle:.2f}°)")
    return result


def crop_spectrum(
    image: CCDData,
    x0_rotated: float,
    y0_rotated: float,
    half_height: int = 40,
    sky_offset: int = 60,
    sky_half_height: int = 20,
) -> tuple[CCDData, CCDData]:
    """
    Crop the rectified image to isolate the spectrum and a sky band.

    After rotation the order 0 is at (x0_rotated, y0_rotated).
    We cut a horizontal strip centred on y0_rotated for the spectrum,
    and a parallel strip offset by sky_offset pixels for the sky.

    Args:
        image (CCDData):    rectified image
        x0_rotated (float): x position of order 0 after rotation
        y0_rotated (float): y position of order 0 after rotation
        half_height (int):  half-height of the spectrum extraction strip (px)
        sky_offset (int):   distance from trace centre to sky strip centre (px)
        sky_half_height (int): half-height of the sky strip (px)

    Returns:
        tuple: (spectrum_strip CCDData, sky_strip CCDData)
               spectrum_strip columns start at order 0 (x0_rotated)
    """
    h, w = image.data.shape
    y_c = int(round(y0_rotated))
    x_start = max(0, int(round(x0_rotated)))

    # spectrum strip
    y_lo = max(0, y_c - half_height)
    y_hi = min(h, y_c + half_height)
    spc_data = image.data[y_lo:y_hi, x_start:].copy()

    # sky strip (above trace)
    sky_y_c = y_c - sky_offset
    sky_lo = max(0, sky_y_c - sky_half_height)
    sky_hi = min(h, sky_y_c + sky_half_height)
    sky_data = image.data[sky_lo:sky_hi, x_start:].copy()

    hdr = image.header.copy()
    hdr['SA100X0']  = int(x_start)
    hdr['SA100Y0']  = y_c
    hdr['SA100YLO'] = y_lo
    hdr['SA100YHI'] = y_hi

    spc_strip = CCDData(spc_data, unit=image.unit, header=hdr)
    sky_strip = CCDData(sky_data, unit=image.unit, header=hdr.copy())

    logging.info(
        f"SA100: spectrum strip y=[{y_lo}:{y_hi}], "
        f"sky strip y=[{sky_lo}:{sky_hi}], "
        f"from x={x_start}"
    )
    return spc_strip, sky_strip


def rotate_point(x: float, y: float, cx: float, cy: float, angle_deg: float) -> tuple[float, float]:
    """
    Rotate a point (x, y) around centre (cx, cy) by angle_deg degrees.

    Used to find where is order 0 after image rectification.

    Args:
        x, y:       original point coordinates
        cx, cy:     rotation centre (usually image centre)
        angle_deg:  rotation angle in degrees (same sign as passed to rectify)

    Returns:
        tuple: (x_rotated, y_rotated)
    """
    rad = math.radians(-angle_deg)
    dx, dy = x - cx, y - cy
    xr = cx + dx * math.cos(rad) - dy * math.sin(rad)
    yr = cy + dx * math.sin(rad) + dy * math.cos(rad)
    return xr, yr


# ---------------------------------------------------------------------------
# Calibration
# ---------------------------------------------------------------------------

def compute_dispersion(distance_mm: float, pixel_size_um: float,
                        ref_wavelength_nm: float = 550.0) -> float:
    """
    Compute the linear dispersion of the SA100/200 in Å/pixel.

    dispersion (Å/px) = (pixel_size_µm x 1000) / (distance_mm x tan(arcsin(λ/d)))

    Uses the grating equation in first order at a reference wavelength:
        sin(θ) = λ / d
        dispersion = pixel_size / (distance x cos(θ))   [same units]

    Args:
        distance_mm (float):     distance from SA100/200 to sensor (mm)
        pixel_size_um (float):   pixel size of the sensor (µm)
        ref_wavelength_nm (float): reference wavelength for linearisation (nm)

    Returns:
        float: dispersion in Å/pixel (positive, increasing toward red)
    """
    d_nm = Config().get_float('sa100', 'grating_period_nm') or 5000.0
    lambda_nm = ref_wavelength_nm
    sin_theta = lambda_nm / d_nm
    cos_theta = math.sqrt(1.0 - sin_theta ** 2)
    # tan(θ) = opposite / adjacent = (n_px × pixel_size) / distance
    # --> dλ/dpx = d × cos(θ) × pixel_size / distance  [all in same unit]
    # pixel_size_um / 1000 --> mm,  result in nm/px × 10 --> Å/px
    dispersion_nm_per_px = (d_nm * cos_theta * (pixel_size_um / 1000.0)) / distance_mm
    dispersion_AA_per_px = dispersion_nm_per_px * 10.0

    logging.info(
        f"SA100: dispersion = {dispersion_AA_per_px:.3f} Å/px "
        f"(distance={distance_mm} mm, pixel={pixel_size_um} µm)"
    )
    return dispersion_AA_per_px


def build_wavelength_axis(
    n_pixels: int,
    dispersion_AA_per_px: float,
    lambda_zero_AA: float = 0.0,
    order_zero_offset_px: int = 0,
) -> np.ndarray:
    """
    Build a linear wavelength array from order 0.

    Args:
        n_pixels (int):              number of pixels along dispersion axis
        dispersion_AA_per_px (float): Å per pixel
        lambda_zero_AA (float):      wavelength at order 0 = 0 Å (by definition)
        order_zero_offset_px (int):  fine-tuning offset in pixels

    Returns:
        np.ndarray: wavelength array in Å
    """
    pixels = np.arange(n_pixels) - order_zero_offset_px
    return lambda_zero_AA + pixels * dispersion_AA_per_px


def calibrate_sa100(
    flux_1d: np.ndarray,
    distance_mm: float,
    pixel_size_um: float,
    order_zero_offset_px: int = 0,
    lambda_min: float = 3500.0,
    lambda_max: float = 7500.0,
) -> Spectrum:
    """
    Build a calibrated Spectrum from a 1D flux array extracted from a SA100 image.

    Args:
        flux_1d (np.ndarray):        extracted flux (ADU), starting at order 0
        distance_mm (float):         SA100-to-sensor distance (mm)
        pixel_size_um (float):       sensor pixel size (µm)
        order_zero_offset_px (int):  fine pixel offset if needed

    Returns:
        Spectrum: calibrated 1D spectrum with wavelength axis in Å
    """
    dispersion = compute_dispersion(distance_mm, pixel_size_um)

    wave = build_wavelength_axis(
        n_pixels=len(flux_1d),
        dispersion_AA_per_px=dispersion,
        order_zero_offset_px=order_zero_offset_px,
    )
    # keep only the requested wavelength range
    mask = (wave > lambda_min) & (wave < lambda_max)
    logging.info(f"SA100: wavelength range [{lambda_min:.0f}–{lambda_max:.0f}] Å, "
                 f"{mask.sum()} pixels kept")
    spc = Spectrum(
        spectral_axis=wave[mask] * u.Angstrom,
        flux=flux_1d[mask] * u.adu,
        )

    cur_time: str = datetime.now().strftime('%m%Y_%H%M%S')
    #spc.write(f"spectrum1D_{cur_time}.fits", overwrite=True)  # format par défaut : tabular-fits

    data = np.column_stack([spc.spectral_axis.value, spc.flux.value])
    np.savetxt(f"spectrum1D_{cur_time}.dat", data,
            header="wavelength(AA)  flux(ADU)",
            comments='# ')

    #tab = Table([spc.spectral_axis.value, spc.flux.value],
     #           names=['wavelength', 'flux'],
      #          meta={'CUNIT1': 'Angstrom', 'CUNIT2': 'ADU'})
    #tab.write(f"spectrum1D_{cur_time}_tbl.dat", format='ascii', overwrite=True)

    return spc


