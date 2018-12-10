#!/usr/bin/env python

import numpy as np
import coastal_tools as ct


def cellWidthVsLatLon():

    km = 1000.0

    params = ct.default_params
    params["dx_min_coastal"] = 20.0*km
    params["trans_width"] = 1200.0*km
    params["trans_start"] = 200.0*km
    params["mesh_type"] = "EC"
    params["dx_max_global"] = 120.0*km
    params["n_longest"] = 20

    params["plot_box"] = np.array([-140,0,30,90]) 

    #print "***Gulf Coast***"
    #params["region_box"] = ct.US_Gulf_Coast
    #params["restrict_box"] = ct.Gulf_restrict
    #cell_width, lon, lat = ct.coastal_refined_mesh(params)

    #print "***East Coast***"
    #params["region_box"] = ct.US_East_Coast
    #params["restrict_box"] = ct.East_Coast_restrict
    #cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    #print "***Caribbean***"
    #params["region_box"] = ct.Caribbean
    #params["restrict_box"] = ct.Caribbean_restrict
    #params["trans_width"] = 400.0*km
    #params["trans_start"] = 300.0*km
    #params["n_longest"] = 50
    #cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)



    print "***Northern Canada***"
    params["region_box"] = {"include":[np.array([-126.5,-59,67,85])],
                            "exclude":[]}
    params["restrict_box"] = ct.Empty
    cell_width, lon, lat = ct.coastal_refined_mesh(params)

    print "***Newfoundland***" 
    Newfoundland = {"include":[np.array([-65.0,-50.0,53.0,60.0])],
                    "exclude":[]}
    params["region_box"] = Newfoundland 
    params["restrict_box"] = ct.Empty
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print "***Greenland***"
    Greenland = {"include":[np.array([-81.5,-12.5,60,85])],
                 "exclude":[np.array([[-87.6,58.7],
                                      [-51.9,56.6],
                                      [-68.9,75.5], 
                                      [-107.0,73.3]]),
                            np.array([[-119.0,74.5],
                                      [-92.7,75.9],
                                      [-52.84,83.25],
                                      [-100.8,84.0]]),
                            np.array([[-101.3,68.5],
                                      [-82.4,72.3],
                                      [-68.7,81.24],
                                      [-117.29,77.75]]),
                            np.array([-25.0,-10.0,62.5,67.5])]}  
    params["region_box"] = Greenland 
    params["restrict_box"] = ct.Empty
    cell_width, lon, lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print "***Bering Straight****"
    params["dx_min_coastal"] = 10.0*km
    params["trans_width"] = 700.0*km
    params["trans_start"] = 100.0*km
    params["lon_wrap"] = True
    params["region_box"] = {"include":[np.array([-172.5,-165,65,66.5])],
                            "exclude":[]} 
    cell_width,lon,lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    print "***Greenland Coast***"
    Greenland = {"include":[np.array([[-50.89,57.89],
                                      [-36.40,59.54],
                                      [-28.31,67.63],
                                      [-59.18,76.02]])],
                 "exclude":[]}
    params["region_box"] = Greenland 
    params["restrict_box"] = ct.Empty
    params["dx_min_coastal"] = 10.0*km
    params["trans_width"] = 200.0*km
    params["trans_start"] = 100.0*km
    params["lon_wrap"] = False  
    cell_width, lon, lat = ct.coastal_refined_mesh(params,cell_width,lon,lat)

    return cell_width / km, lon, lat
