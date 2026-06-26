import os
import speclite.filters
from FilterMags.filtermag.__init__ import DATADIR
import FilterMags.filtermag as filtermag
import math
import unittest

def test_read_spectrum():
    """
    Unit test for ``read_spectrum`` standard input.
    """
    spectrumfile0 = os.path.join(DATADIR, 'spec-5164-56067-0822.fits')
    wave0, flux0 = filtermag.read_spectrum(spectrumfile0)
    assert len(wave0) == len(flux0)

    spectrumfile1 = os.path.join(DATADIR, 'spec1-1678-53433-0425.fits')
    wave1, flux1 = filtermag.read_spectrum(spectrumfile1)
    assert len(wave1) == len(flux1)

    spectrumfile2 = os.path.join(DATADIR, 'spec2-1678-53433-0429.fits')
    wave2, flux2 = filtermag.read_spectrum(spectrumfile2)
    assert len(wave2) == len(flux2)

test_read_spectrum()

def test_filter_mag():
    """
    End-to-end test for ``filter_mag`` standard input.
    """
    spectrumfile0 = os.path.join(DATADIR, 'spec-5164-56067-0822.fits')
    wave0, flux0 = filtermag.read_spectrum(spectrumfile0)
    rmag0 = filtermag.filter_mag(wave0, flux0, 'sdss_r')
    rband = speclite.filters.load_filter('sdss2010-r')
    astropy_rmag0 = rband.get_ab_magnitude(flux0, wave0)
    assert math.isclose(rmag0, astropy_rmag0, rel_tol=0.01)

    spectrumfile1 = os.path.join(DATADIR, 'spec1-1678-53433-0425.fits')
    wave1, flux1 = filtermag.read_spectrum(spectrumfile1)
    rmag1 = filtermag.filter_mag(wave1, flux1, 'sdss_r')
    astropy_rmag1 = rband.get_ab_magnitude(flux1, wave1)
    assert math.isclose(rmag1, astropy_rmag1, rel_tol=0.01)

    spectrumfile2 = os.path.join(DATADIR, 'spec2-1678-53433-0429.fits')
    wave2, flux2 = filtermag.read_spectrum(spectrumfile2)
    iband = speclite.filters.load_filter('sdss2010-i')
    imag2 = filtermag.filter_mag(wave2, flux2, 'sdss_i')
    astropy_imag2 = iband.get_ab_magnitude(flux2, wave2)
    assert math.isclose(imag2, astropy_imag2, rel_tol=0.01)
    
test_filter_mag()

class TestWithErrors(unittest.TestCase):
    """
    Class for unit tests of edge cases meant to throw errors.
    """
    def __init__(self):
        pass

    def test_out_of_range(self):
        """
        Unit test for when the spectrum wavelength is incompatible with the transmission curve wavelength
        range.
        """
        spectrumfile = os.path.join(DATADIR, 'spec-5164-56067-0822.fits')
        wave, flux = filtermag.read_spectrum(spectrumfile)
        with self.assertRaises(ValueError):
            filtermag.filter_mag(wave, flux, 'sdss_z')
        pass

out_of_range = TestWithErrors().test_out_of_range()
