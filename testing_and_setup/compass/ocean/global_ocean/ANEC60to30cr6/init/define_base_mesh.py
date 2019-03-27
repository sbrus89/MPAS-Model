#!/usr/bin/env python

import numpy as np
import coastal_tools as ct


def cellWidthVsLatLon():

    km = 1000.0

    params = ct.default_params
    params["dx_min_coastal"] = 20.0*km
    params["trans_width"] = 1200.0*km
    params["trans_start"] = 200.0*km
    params["mesh_type"] = "QU"
    params["dx_max_global"] = 120.0*km
    params["nc_file"] = './ETOPO2v2c_f4_151106.nc'
    params["nc_vars"] = ["x","y","z"]

    params["plot_box"] = np.array([-180,180,-80,-40]) 


    params["region_box"] = {"include":[np.array([-180,180,-80,-60])],
                            "exclude":[]}
    params["restrict_box"] = ct.Empty
    params["z_contour"] = -900.0
    cell_width, lon, lat = ct.coastal_refined_mesh(params)


    return cell_width / km, lon, lat
