''' 

Ps:saturated soil water potential; MPa
    b: exponent on the water retention curve
    n: porosity
    Ksat: saturated soil hydraulic conductivity; m/d; 
    
http://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/survey/office/ssr10/tr/?cid=nrcs144p2_074846
Parameters for loamy_sandy derived from measurements at field sites in Morongo Valley

'''
soil_dict = {'sand':        {'Ps': -0.34*10**(-3), 'b':4.05, 'n':0.35, 'Ksat':2.0, 'sh':0.08, 'sfc':0.35, 'sst':0.33},
             'loamy_sand':  {'Ps': -0.17*10**(-3), 'b':4.38, 'n':0.42, 'Ksat':1.0, 'sh':0.08, 'sfc':0.52, 'sst':0.31},
             'sandy_loam':  {'Ps': -0.70*10**(-3), 'b':4.90, 'n':0.43, 'Ksat':0.80, 'sh':0.14, 'sfc':0.56, 'sst':0.46},
             'loam':        {'Ps': -1.43*10**(-3), 'b':5.39, 'n':0.45, 'Ksat':0.20, 'sh':0.19, 'sfc':0.65, 'sst':0.57},
             'clay':        {'Ps': -1.82*10**(-3), 'b':11.4, 'n':0.50, 'Ksat':0.04, 'sh':0.47, 'sfc':1.0, 'sst':0.78}}

