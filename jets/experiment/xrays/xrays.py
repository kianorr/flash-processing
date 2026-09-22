import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from skimage import exposure, transform, data
import os
import glob
from scipy.ndimage import gaussian_filter, median_filter, uniform_filter
import h5py
from numpy.lib.stride_tricks import sliding_window_view
# from flashtools.utils import get_closest
from matplotlib.animation import FuncAnimation
import pickle
import pandas as pd
from scipy.signal import find_peaks
from scipy.interpolate import PchipInterpolator
from scipy.optimize import curve_fit
import warnings


def load_tiff(
    shot,
    x_resolution=314572800,
    y_resolution=262144,
    magnification=1.5,
    foreshortening=None,
    diagnostic='xrfc3',
):
    dpi = x_resolution / y_resolution
    inches_to_mm = 25.4
    pixel_size_mm = inches_to_mm / dpi
    if magnification:
        pixel_size_mm /= magnification
    if foreshortening:
        pixel_size_mm /= foreshortening
    path = "/Users/johan/Documents/research/data_analysis/XRFC"
    file = f'XRFC-s{shot}_{diagnostic}_-1.tif'
    image = np.array(Image.open(os.path.join(path, file)))
    image = transform.rotate(image, angle=90, resize=True)[:, :, 0]
    xaxis = np.arange(np.shape(image)[1]) * pixel_size_mm
    yaxis = np.arange(np.shape(image)[0]) * pixel_size_mm

    return xaxis, yaxis, image


def load_pds(
    shot,
    magnification=1.5,
    foreshortening=None,
):
    path = "/Users/johan/Documents/research/data_analysis/XRFC"
    wfile = f"XRFC3-xrf3t1_{shot}_swedge.h5"
    with h5py.File(os.path.join(path, wfile), "r") as f:
        wedge = np.array(f['pds_image'])

    with h5py.File(os.path.join(path, f"XRFC3-xrf3t1_{shot}.h5"), "r") as f:
        data = np.array(f['pds_image'])
        dim0 = np.array(f['fakeDim0'])
        dim1 = np.array(f['fakeDim1'])
    full_image = convert_wedge_data(data, np.mean(wedge, axis=0), visualize=False)
    conversion = 1.
    if magnification:
        conversion /= magnification
    if foreshortening:
        conversion /= foreshortening
    xaxis = dim0 * conversion * 1e-3
    yaxis = dim1 * conversion * 1e-3

    return xaxis, yaxis, full_image.T


def load_ccd(
    shot,
    pixel_size_um=9,
    magnification=1.5,
    binning=2,
    foreshortening=None,
):

    if binning:
        pixel_size_um *= binning
    if magnification:
        pixel_size_um /= magnification
    if foreshortening:
        pixel_size_um /= foreshortening

    pixel_size_mm = pixel_size_um * 1e-3
    path = "/Users/johan/Documents/research/data_analysis/XRFC"
    file = f"XRFC3CCD-xrfc3_t1_{shot}.h5"
    with h5py.File(os.path.join(path, file), "r") as f:
        data = np.array(f['Streak_array'])
        dim0 = np.array(f['fakeDim0'])
        dim1 = np.array(f['fakeDim1'])
        dim2 = np.array(f['fakeDim2'])

    ccd_image = transform.rotate(data[0], angle=-90, resize=True)
    yaxis = np.arange(np.shape(ccd_image)[0]) * pixel_size_mm
    xaxis = np.arange(np.shape(ccd_image)[1]) * pixel_size_mm

    return xaxis, yaxis, ccd_image[::-1]


def split_into_frames(
    xaxis, yaxis, image,
    target_sep = 6.4,
    target_thickness=0.1,
    y_center=19.5,
    x_centers=np.array([2.5, 8.5, 14.7, 20.7]),
    x_radius=1.5,
):
    yaxis_min = y_center - (target_sep / 2) - target_thickness
    yaxis_max = yaxis_min + target_thickness * 2 + target_sep
    yrange = (yaxis > yaxis_min) & (yaxis < yaxis_max)
    new_yaxis = yaxis[yrange]
    y_radius = (yaxis_max - yaxis_min) / 2
    new_yaxis = (new_yaxis - (np.min(new_yaxis) + np.max(new_yaxis)) / 2)


    new_xaxis = (xaxis - np.max(xaxis) / 2)
    xrange = (new_xaxis > -x_radius) & (new_xaxis < x_radius)
    new_xaxis = new_xaxis[xrange]


    center_inds = get_closest(xaxis, x_centers)
    xaxis_center_ind = get_closest(xaxis, x_centers[0])

    xaxis_max = x_centers + x_radius
    xaxis_max_ind = get_closest(xaxis, xaxis_max)
    difference = abs(center_inds[0] - xaxis_max_ind[0])

    xaxis_slice = [
        slice(center_inds[i] - difference, center_inds[i] + difference)
        for i in range(len(x_centers))
    ]

    frames = np.zeros((len(x_centers), len(new_yaxis), difference * 2))
    for i, center in enumerate(center_inds):
        frames[i] = image[yrange, xaxis_slice[i]]

    return new_xaxis, new_yaxis, frames


def window_average_frames(xaxis, yaxis, frames, xaxis_window=1, yaxis_window=1):
    frames_avg = sliding_window_view(frames, window_shape=xaxis_window, axis=2).mean(axis=-1)
    xaxis_avg = sliding_window_view(xaxis, window_shape=xaxis_window).mean(axis=-1)
    
    frames_avg = sliding_window_view(frames_avg, window_shape=yaxis_window, axis=1).mean(axis=-1)
    yaxis_avg = sliding_window_view(yaxis, window_shape=yaxis_window).mean(axis=-1)

    return xaxis_avg, yaxis_avg, frames_avg


def convert_wedge_data(
    data,
    wedge,
    visualize=True,
    nsteps=7,
    edge_buffer=0.3,
    nwedge=21,
    atten_start=0.0,
    atten_step=0.15,
    fit_type='sqrt((a*x+b)**2 + c**2) - a*x + d',
    fit_start=(100, -1500, 400, 1500)
):
    """
    Converts PDS-scanned images to signal level.

    Parameters
    ----------
    data : 2D array
        The image to be converted.
    wedge : 1D or 2D array
        The associated wedge image (will be averaged vertically).
    visualize : bool, optional
        Show conversion plots. Default is True.
    nsteps : int, optional
        Number of wedge steps from which to calculate distances. Default is 7.
    edge_buffer : float, optional
        Fraction to shrink between the identified valleys. Default is 0.3.
    nwedge : int, optional
        Total number of steps. Default is 21.
    atten_start : float, optional
        Attenuation of the brightest wedge. Default is 0.
    atten_step : float, optional
        Increment of subsequent wedges. Default is 0.15.
    fit_type : str, optional
        Functional form for curve fitting. Set to 'none' to disable fitting.
    fit_start : tuple, optional
        Starting parameters for the fit function.

    Returns
    -------
    converted : 2D array
        Converted signal level image.
    """

    # --- ensure double precision ---
    data = np.asarray(data, dtype=float)
    wedge = np.asarray(wedge, dtype=float)

    # --- vertically average the wedge ---
    if wedge.ndim > 1 and wedge.shape[0] > 1:
        wedge = wedge.mean(axis=0)

    # --- find valleys (minima between wedge steps) ---
    minprom = np.max(wedge) / 20
    peaks, _ = find_peaks(-wedge, prominence=minprom)

    if len(peaks) < nsteps + 1:
        raise ValueError("Not enough wedge minima found for the requested nsteps.")

    x0 = peaks[0]
    dx = (peaks[nsteps] - peaks[0]) / nsteps

    # shift if there's space to the left of first peak
    if x0 - dx * (1 - edge_buffer) > 0:
        x0 -= dx

    # --- compute signal levels between valleys ---
    signal = np.zeros(nwedge)
    isused = np.zeros_like(wedge, dtype=bool)
    for i in range(1, nwedge+1):
        xmin = x0 + (i - 1 + edge_buffer) * dx if i > 0 else x0 + edge_buffer * dx
        xmax = x0 + (i - edge_buffer) * dx
        xrange = np.arange(int(np.floor(xmin)), int(np.ceil(xmax)))
        xrange = xrange[(xrange >= 0) & (xrange < len(wedge))]

        if len(xrange) == 0:
            continue

        signal[i-1] = np.mean(wedge[xrange])

        # enforce strictly decreasing
        # if i > 0:
        #     signal[i-1] = min(signal[i-1], 0.99 * signal[i - 2])
        # print(signal)
        isused[xrange] = True
    # --- attenuation and transmission ---
    atten = atten_start + np.arange(nwedge) * atten_step
    trans = 10 ** (-atten)

    # --- define fit function ---
    def fit_func(x, a, b, c, d):
        return np.sqrt((a * x + b) ** 2 + c**2) - a * x + d

    # --- fit or interpolate ---
    if fit_type.lower() == "none":
        f = lambda x: np.interp(x, atten, signal)
    else:
        try:
            popt, _ = curve_fit(fit_func, atten, signal, p0=fit_start, maxfev=10000)
            f = lambda x: fit_func(np.array(x), *popt)
        except Exception as e:
            warnings.warn(f'Poor fit ({e}). Reverting to raw points.')
            f = lambda x: np.interp(x, atten, signal)
            fit_type = "none"

    signal_fit = f(atten)

    # --- ensure no extrapolation below last signal ---
    data_floor = np.maximum(data, signal_fit[-1])
    # data_floor = data
    # --- convert data via interpolation ---
    # if signal_fit[0] > signal_fit[-1]:
    #     signal_fit = signal_fit[::-1]
    #     atten = atten[::-1]
        # trans=trans[::-1]
    interp = PchipInterpolator(signal_fit[::-1], atten[::-1])
    # interp = interp1d(signal_fit, atten)
    converted = interp(data_floor)
    converted = 10 ** (-converted)

    # --- visualization ---
    if visualize:
        fig, axs = plt.subplots(2, 3, figsize=(14, 8))
        fig.suptitle("Convert Wedge Data", fontsize=16, fontweight='bold')

        # Input data
        im = axs[0, 0].imshow(data, cmap='gray')
        axs[0, 0].set_title("Input Data")
        axs[0, 0].axis("off")
        fig.colorbar(im, ax=axs[0, 0], label="Signal")

        # Wedge data
        im = axs[1, 0].imshow(np.atleast_2d(wedge), cmap='gray', aspect='auto')
        axs[1, 0].set_title("Wedge Data")
        axs[1, 0].axis("off")
        fig.colorbar(im, ax=axs[1, 0], label="Signal")

        # Wedge extraction
        axs[0, 1].plot(wedge, 'k:', label='Wedge Lineout')
        used_wedge = np.copy(wedge)
        used_wedge[~isused] = np.nan
        axs[0, 1].plot(used_wedge, 'k-', lw=2, label='Used Data')
        axs[0, 1].plot(x0 + (np.arange(nwedge) + 0.5) * dx, signal, 'rx', lw=2, label='Estimated Signal')

        if fit_type.lower() != "none":
            fx = (((np.arange(len(wedge)) - x0) / dx) - 0.5) * atten_step + atten_start
            axs[0, 1].plot(f(fx), 'b--', lw=1, label=f'Fitted Signal ({fit_type})')

        axs[0, 1].set_title("Wedge Extraction")
        axs[0, 1].set_xlabel("Position (pixel)")
        axs[0, 1].set_ylabel("Signal")
        axs[0, 1].legend()
        axs[0, 1].grid(True)

        # Wedge conversion
        axs[1, 1].plot(signal_fit, atten, '-o')
        axs[1, 1].set_ylabel("Attenuation", color='C0')

        x_vals = np.linspace(signal_fit[-1], np.max(signal_fit), 200)
        y_vals = PchipInterpolator(signal_fit[::-1], trans[::-1])(x_vals)
        y_vals = np.maximum(0, y_vals)
        ax2 = axs[1, 1].twinx()
        ax2.plot(signal_fit, trans, 'o', color='orange')
        ax2.plot(x_vals, y_vals, '-', color='orange')
        ax2.set_ylabel("Transmission", color='orange')
        axs[1, 1].set_xlabel("Signal")
        axs[1, 1].set_title("Wedge Conversion")

        # Converted data
        im = axs[0, 2].imshow(converted, cmap='gray', norm="log")
        axs[0, 2].set_title("Results")
        axs[0, 2].axis("off")
        fig.colorbar(im, ax=axs[0, 2], label="Converted Data")

        plt.tight_layout()
        plt.show()

    return converted
