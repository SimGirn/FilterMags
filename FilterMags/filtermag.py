# outside packages
import os
from astropy.io import ascii
import numpy as np

# dictionary containing keywords for all filter transmission files
current_dir = os.getcwd()
filter_dict = {
    "sdss_u": os.path.join(current_dir, 'TransmissionCurves/SLOAN_SDSS.u.dat'),
    "sdss_g": os.path.join(current_dir, 'TransmissionCurves/SLOAN_SDSS.g.dat'),
    "sdss_r": os.path.join(current_dir, 'TransmissionCurves/SLOAN_SDSS.r.dat'),
    "sdss_i": os.path.join(current_dir, 'TransmissionCurves/SLOAN_SDSS.i.dat'),
    "sdss_z": os.path.join(current_dir, 'TransmissionCurves/SLOAN_SDSS.z.dat'),
}

# read the transmission curve
class TransmissionCurve:
    """
    Class that defines transmission curve objects. Each transmission curve has a filter name and arrays of 
    wavelengths and transmission (%) based on the transmission curve files in the repo. The dictionary 
    ``filter_dict`` contains the file names for each filter.
    """
    def __init__(self, filtername):
        """Constructor method
        """
        filename = filter_dict[filtername]
        transmission_data = ascii.read(filename)
        self.wavelength = transmission_data['col1'].data
        self.transmission = transmission_data['col2'].data
        pass

# calculate filter magnitude
def filter_mag(wave_spectra, flux, filter):
    """Filter magnitude calculator

    Calculate the filter magnitude for a spectrum provided by the user in a given filter.
    
    Args: 
        wave_spectra (array): the wavelengths of the spectrum (in Angstroms)
        flux (array): the associated fluxes for the wavelengths (in erg s-1 cm-2 A-1)
        filter (string): the filter that magnitude will be calculated in, identified using a key in ``filter_dict``
    
    Returns: 
        float: filter magnitude for the spectrum
    """
    # read in the transmission curve for the given filter
    trans_curve = TransmissionCurve(filter)
    wave_trans = trans_curve.wavelength
    t = trans_curve.transmission

    if max(wave_spectra) < max(wave_trans) or min(wave_spectra) > min(wave_trans):
        raise ValueError('At least part of filter out of range of spectrum')
    
    # three integrals that are needed for calculation
    int_1 = np.trapezoid(t / wave_trans, wave_trans)
    int_2 = np.trapezoid(t * wave_trans, wave_trans)
    t_spectra = np.interp(wave_spectra, wave_trans, t, left=0, right=0) #interpolation for third integral
    int_3 = np.trapezoid((flux*t_spectra)/wave_spectra, wave_spectra)

    # magnitude calculation
    flam = int_3 / int_1 #f_lambda
    pivot = int_2 / int_1 #pivot wavelength squared
    c = 3e18 #angstrom/s
    mag = -2.5 * np.log10(flam) - 2.5 * np.log10(pivot / c) - 48.6

    return mag
