#!/usr/bin/python

from astropy.io import fits
import numpy as np
import healpy as hp
import sys
import math
import matplotlib.pyplot as plt
import scipy.io
from scipy.interpolate import griddata


class ImageFromFits:

    def __init__(self, signal_arr, ra_axis=None, dec_axis=None,
                 ra_range=None, dec_range=None):
        n_dec_vals, n_ra_vals = np.shape(signal_arr)
        if ra_axis is not None:
            ra_axis = list(ra_axis)
            if len(ra_axis) != n_ra_vals:
                print('ERROR: Number of axis elements does not match data axis. Exiting.')
                sys.exit(1)
        if dec_axis is not None:
            dec_axis = list(dec_axis)
            if len(dec_axis) != n_dec_vals:
                print('ERROR: Number of axis elements does not match data axis. Exiting.')
                sys.exit(1)
        if ra_axis is None and ra_range is not None:
            if len(ra_range) != 2:
                print('ERROR: parameters ra_range and dec_range must have the form [min, max]. Exiting.')
                sys.exit(1)
            ra_axis = list(np.linspace(ra_range[0], ra_range[1], n_ra_vals))
        if dec_axis is None and dec_range is not None:
            if len(dec_range) != 2:
                print('ERROR: parameters ra_range and dec_range must have the form [min, max]. Exiting.')
                sys.exit(1)
            dec_axis = list(
                np.linspace(dec_range[0], dec_range[1], n_dec_vals)
            )
            print('uh oh')
        self.signal_arr = signal_arr
        self.ra_axis = ra_axis
        self.dec_axis = dec_axis

    def limit_data_range(self, ra_range=None, dec_range=None):
        if ra_range is not None:
            use_ra_inds = [i for i in range(len(self.ra_axis))
                           if ra_range[0] < self.ra_axis[i] < ra_range[1]
                           ]
        else:
            use_ra_inds = range(len(self.ra_axis))
            print('use_ra_inds:' + use_ra_inds)
        if dec_range is not None:
            use_dec_inds = [i for i in range(len(self.dec_axis))
                            if dec_range[0] < self.dec_axis[i] < dec_range[1]
                            ]
            print('whoops')
        else:
            use_dec_inds = range(len(self.dec_axis))
            print('use_dec_inds:' + use_dec_inds)
        self.signal_arr = self.signal_arr[
            use_dec_inds[0]:use_dec_inds[-1]+1,
            use_ra_inds[0]:use_ra_inds[-1]+1
        ]
        self.ra_axis = list(
            np.linspace(ra_range[0], ra_range[1], len(use_ra_inds))
        )
        self.dec_axis = list(
            np.linspace(dec_range[0], dec_range[1], len(use_dec_inds))
        )


def load_image(data_filename):

    contents = fits.open(data_filename)
    use_hdu = 0
    data = contents[use_hdu].data
    header = contents[use_hdu].header

    if 'CD1_1' in header.keys() and 'CD2_2' in header.keys():  # FHD convention
        cdelt1 = header['CD1_1']
        cdelt2 = header['CD2_2']
        print ('WARNING: Ignoring curved sky effects.')
    elif 'CDELT1' in header.keys() and 'CDELT2' in header.keys():
        cdelt1 = header['CDELT1']
        cdelt2 = header['CDELT2']
    else:
        print('ERROR: Header format not recognized.')
        sys.exit(1)

    ra_axis = [
        header['crval1'] +
        cdelt1*(i-header['crpix1'])
        for i in range(header['naxis1'])
        ]
    dec_axis = [
        header['crval2'] +
        cdelt2*(i-header['crpix2'])
        for i in range(header['naxis2'])
        ]

    fits_image = ImageFromFits(data, ra_axis=ra_axis, dec_axis=dec_axis)
    return fits_image


def load_gaussian_source_model_as_image(
    catalog_path, source_ind=0, resolution=.01, ra_range=None, dec_range=None,
    reference_image=None
):
    # THE NORMALIZATION IS WRONG!!! DON'T USE THIS FUNCTION!

    # If reference image is supplied, use the same ra and dec locations
    if reference_image is not None:
        print('ref image is not none')
        ra_axis = reference_image.ra_axis
        dec_axis = reference_image.dec_axis
        grid_ra = np.tile(np.array(ra_axis), (len(dec_axis), 1))
        grid_dec = np.tile(np.array([dec_axis]).T, (1, len(ra_axis)))
        ra_range = [min(ra_axis), max(ra_axis)]
        dec_range = [min(dec_axis), max(dec_axis)]
    else:
        if ra_range is None:
            ra_range = [50, 51.25]
        if dec_range is None:
            dec_range = [-37.8, -36.7]

        grid_dec, grid_ra = np.mgrid[
            dec_range[0]:dec_range[1]:resolution,
            ra_range[0]:ra_range[1]:resolution
            ]

    plot_signal = np.zeros_like(grid_dec)

    source = scipy.io.readsav(catalog_path)['catalog'][source_ind]
    source_ra = source['ra']
    source_dec = source['dec']
    components = source['extend']
    if len(components) == 0:
        print('WARNING: Source is not extended.')

    total_flux = 0.
    for comp in components:
        comp_ra = comp['ra']
        comp_dec = comp['dec']
        comp_flux = comp['flux']['I'][0]
        total_flux += comp_flux
        print('total flux is' + total_flux)
        comp_size_x = comp['shape']['x'][0]/(7200.*np.sqrt(2*np.log(2.)))
        print(comp_size_x)
        comp_size_y = comp['shape']['y'][0]/(7200.*np.sqrt(2*np.log(2.)))
        comp_size_angle = comp['shape']['angle'][0]

        if comp_size_x == 0:
            comp_size_x = resolution
        if comp_size_y == 0:
            comp_size_y = resolution

        for i in range(np.shape(grid_dec)[0]):
            for j in range(np.shape(grid_dec)[1]):
                pixel_val = (
                    comp_flux*np.pi/(2.*(180.)**2.*comp_size_x*comp_size_y)
                    * np.exp(-(grid_ra[i, j]-comp_ra)**2./(2*comp_size_x**2.))
                    * np.exp(-(grid_dec[i, j]-comp_dec)**2./(2*comp_size_y**2.))
                )
                plot_signal[i, j] += pixel_val

    image = ImageFromFits(plot_signal, ra_range=ra_range, dec_range=dec_range)
    return image


def difference_images(image1, image2):

    if image1.ra_axis == image2.ra_axis and image1.dec_axis == image2.dec_axis:
        data_diff = ImageFromFits(
            np.subtract(image1.signal_arr, image2.signal_arr),
            ra_axis=image1.ra_axis, dec_axis=image1.dec_axis
            )
    else:
        print ('WARNING: Image axes do not match. Interpolating image2 to image1 axes.')
        image2_signal_array_interp = griddata(  # This doesn't work
            (image2.ra_axis, image2.dec_axis),
            image2.signal_arr,
            (image1.ra_axis, image1.dec_axis)
            )
        data_diff = ImageFromFits(
            np.subtract(image1.signal_arr, image2_signal_array_interp),
            ra_axis=image1.ra_axis, dec_axis=image1.dec_axis
            )
    return data_diff


def plot_fits_image(
    fits_image, color_scale, output_path, title='', ra_range=None, dec_range=None, log=False,
    colorbar_label='Flux Density (Jy/sr)', plot_grid=True,
    xlabel='RA (deg.)', ylabel='Dec. (deg.)'
):

    ra_range = [25,95]
    dec_range = [-65,5]
    colorbar_range = color_scale
    save_filename = output_path
    if ra_range is not None or dec_range is not None:
        fits_image.limit_data_range(ra_range=ra_range, dec_range=dec_range)

    print('Minimum image value {}'.format(np.min(fits_image.signal_arr)))
    print('Maximum image value {}'.format(np.max(fits_image.signal_arr)))

    fig, ax = plt.subplots()
    plt.imshow(
        fits_image.signal_arr, origin='lower', interpolation='none',
        cmap='Greys_r',
        extent=[
            fits_image.ra_axis[0], fits_image.ra_axis[-1],
            fits_image.dec_axis[0], fits_image.dec_axis[-1]
            ],
        vmin=colorbar_range[0], vmax=colorbar_range[1], aspect='auto'
    )
    print(fits_image.ra_axis[0])
    print(fits_image.ra_axis[1])
    plt.axis('equal')
    #ax.set_facecolor('gray')  # make plot background gray
    ax.set_facecolor('black')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if plot_grid:
        plt.grid(which='both', zorder=10, lw=0.5)
    cbar = plt.colorbar()
    # Label colorbar:
    cbar.ax.set_ylabel(colorbar_label, rotation=270, labelpad=15)
    if save_filename == '':
        plt.show()
    else:
        print('Saving figure to {}'.format(save_filename))
        plt.savefig(save_filename, format='png', dpi=500)


if __name__ == '__main__':

    #orig = load_image('/Users/ruby/Downloads/1131478056_uniform_Dirty_XX_decon.fits')
    #new = load_image('/Users/ruby/Downloads/1131478056_uniform_Dirty_XX.fits')
    #diff = difference_images(orig, new)
    print('testing')
    prefix = '/Users/home/Documents/_Files/Dara/School/Graduate/RadCos/FHD_Pyuvsim_comparison/Results/'
    filename = 'ref_1.1_uniform_uniform_Dirty_XX'
    data = load_image(prefix + filename + '.fits')
    output_path = prefix + filename + '_adjusted.png'
    print('testing')
    color_scale = [-.3,1.1]
    plot_fits_image(data, color_scale, output_path)
