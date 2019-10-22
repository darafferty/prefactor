#! /usr/bin/env python
"""
Script to separate a sky model into outlier and field parts
"""
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import casacore.tables as pt
import lsmtool


def main(skymodel, ms_in, fwhm_deg, outroot, scale_factor=1.0, pointing_ra=None,
         pointing_dec=None):
    """
    Separate a makesourcedb sky model into outlier and field parts

    Parameters
    ----------
    skymodel : str
        Filename of the input makesourcedb sky model
    ms_in : str
        Filename of the input MS file
    fwhm_deg : str
        Size of FWHM of primary beam in degrees as "FWHM_RA FWHM_Dec"
    outroot : str
        Root for output sky models: outroot.outlier and outroot.field)
    scale_factor : float
        Scaling to use to determine field region: FWHM * scale_factor
    pointing_ra : float
        RA of pointing in degrees (if None, the phase center is used)
    pointing_dec : float
        Dec of pointing in degrees (if None, the phase center is used)
    """
    ms_in = ms_in.strip('[]')
    fwhm_ra_deg = float(fwhm_deg.split(" ")[0])
    fwhm_dec_deg = float(fwhm_deg.split(" ")[1])
    print('Using width in RA of {0} deg and in Dec of {1} deg'.format(fwhm_ra_deg*2.0, fwhm_dec_deg*2.0))
    scale_factor = float(scale_factor)

    # Get pointing info
    obs = pt.table(ms_in+'::FIELD', ack=False)
    if pointing_ra is None:
        pointing_ra = np.degrees(float(obs.col('REFERENCE_DIR')[0][0][0]))
    if pointing_ra < 0.:
        pointing_ra = 360.0 + (pointing_ra)
    if pointing_dec is None:
        pointing_dec = np.degrees(float(obs.col('REFERENCE_DIR')[0][0][1]))
    obs.close()

    # Filter (in image coordinates) on rectangle defined by FWHM
    s = lsmtool.load(skymodel)
    s_outlier = s.copy()
    s_field = s.copy()
    x, y, refRA, refDec = s._getXY()
    pointing_x, pointing_y = lsmtool.operations_lib.radec2xy([pointing_ra], [pointing_dec],
                                                             refRA=refRA, refDec=refDec)
    crdelt = 0.066667  # deg/pix used by lsmtool.operations_lib.radec2xy()
    min_x = pointing_x[0] - fwhm_ra_deg / crdelt * scale_factor
    max_x = pointing_x[0] + fwhm_ra_deg / crdelt * scale_factor
    min_y = pointing_y[0] - fwhm_dec_deg / crdelt * scale_factor
    max_y = pointing_y[0] + fwhm_dec_deg / crdelt * scale_factor
    field_ind = (x > min_x) & (y > min_y) & (x < max_x) & (y < max_y)

    s_outlier.remove(field_ind)
    s_outlier.write(outroot+'.outlier', clobber=True)
    s_field.select(field_ind)
    s_field.write(outroot+'.field', clobber=True)


if __name__ == '__main__':
    descriptiontext = "Make a makesourcedb sky model from WSClean fits model images.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('fits_models', help='Model images')
    parser.add_argument('msfile', help='MS file')
    parser.add_argument('skymodel', help='Filename of output sky model')
    parser.add_argument('fits_masks', help='Mask images')

    args = parser.parse_args()
    main(args.fits_models, args.msfile, args.skymodel, args.fits_masks)
