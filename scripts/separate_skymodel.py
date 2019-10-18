#! /usr/bin/env python
"""
Script to separate a sky model into outlier and field parts
"""
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import casacore.tables as pt
import lsmtool
from lofarpipe.support.data_map import DataMap


def main(skymodel, ms_input, outroot, scale_factor=1.25):
    """
    Separate a makesourcedb sky model into outlier and field parts

    Parameters
    ----------
    skymodel : str
        Filename of the input makesourcedb sky model
    ms_input : lsit
        List of MS files
    outroot : str
        Root for output sky models: outroot.outlier and outroot.field)
    scale_factor : float
        Scaling to use to determine field region: FWHM * scale_factor
    """
    if type(ms_input) is str:
        if ms_input.startswith('[') and ms_input.endswith(']'):
            ms_list = [f.strip(' \'\"') for f in ms_input.strip('[]').split(',')]
        else:
            map_in = DataMap.load(ms_input)
            map_in.iterator = DataMap.SkipIterator
            ms_list = []
            for fname in map_in:
                if fname.startswith('[') and fname.endswith(']'):
                    for f in fname.strip('[]').split(','):
                        ms_list.append(f.strip(' \'\"'))
                else:
                    ms_list.append(fname.strip(' \'\"'))
    elif type(ms_input) is list:
        ms_list = [str(f).strip(' \'\"') for f in ms_input]
    else:
        raise TypeError('separate_skymodel: type of "ms_input" unknown!')
    scale_factor = int(scale_factor)

    # Find size of primary beam at central frequency
    msfreqs = []
    for ms in ms_list:
        # group all MSs by frequency
        sw = pt.table(ms+'::SPECTRAL_WINDOW', ack=False)
        msfreqs.append(int(sw.col('REF_FREQUENCY')[0]))
    mid_freq = np.mean(msfreqs)
    ant = pt.table(ms_list[0]+'::ANTENNA', ack=False)
    diam = float(ant.col('DISH_DIAMETER')[0])
    ant.close()
    el_values = pt.taql("SELECT mscal.azel1()[1] AS el from "
                        + ms_list[0] + " limit ::10000").getcol("el")
    mean_el_rad = np.mean(el_values)
    sec_el = 1.0 / np.sin(mean_el_rad)
    fwhm_deg_ra = 1.1 * ((3.0e8 / mid_freq) / diam) * 180. / np.pi * sec_el
    fwhm_deg_dec = 1.0 / np.sin(mean_el_rad)

    # Get pointing info
    obs = pt.table(ms_list[0]+'::FIELD', ack=False)
    pointing_ra = np.degrees(float(obs.col('REFERENCE_DIR')[0][0][0]))
    if pointing_ra < 0.:
        pointing_ra = 360.0 + (pointing_ra)
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
    min_x = pointing_x[0] - fwhm_deg_ra / crdelt * scale_factor
    max_x = pointing_x[0] + fwhm_deg_ra / crdelt * scale_factor
    min_y = pointing_y[0] - fwhm_deg_dec / crdelt * scale_factor
    max_y = pointing_y[0] + fwhm_deg_dec / crdelt * scale_factor
    field_ind = (x > min_x) & (y > min_y) & (x < max_x) & (y < max_y)

    s_outlier.remove(field_ind)
    s_outlier.write(outroot+'.outlier')
    s_field.select(field_ind)
    s_field.write(outroot+'.field')


if __name__ == '__main__':
    descriptiontext = "Make a makesourcedb sky model from WSClean fits model images.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('fits_models', help='Model images')
    parser.add_argument('msfile', help='MS file')
    parser.add_argument('skymodel', help='Filename of output sky model')
    parser.add_argument('fits_masks', help='Mask images')

    args = parser.parse_args()
    main(args.fits_models, args.msfile, args.skymodel, args.fits_masks)
