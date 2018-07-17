import xlrd
import numpy as np
from soil_plant_model import Canopy, Stem, Soil_root, Whole_plant, Soil

def import_traits_data(filepath=None, species=None, sp_coln=None):
    if filepath == None: filepath = '../hydraulic_traits.xls'
    if species == None: species = ['JUNI', 'PINE']
    if sp_coln == None: sp_coln = [2,4]
    book = xlrd.open_workbook(filepath)
    sheet = book.sheet_by_name('parameters')
    keys = np.asarray(filter(None, sheet.col_values(0)), dtype='str')
    canopy_keys = ['A_canopy', 'Gs_leaf', 'Amax', 'rho', 'lamABA', 'm']
    stem_keys = ['L_stem', 'A_stem', 'Ksat_stem', 'a_stem', 'P50_stem']
    root_keys = ['L_root', 'A_root', 'd_root']
    chap_dict = {}
    for sp, spc in zip(species, sp_coln):
        chap_dict[sp] = {}
        for part_keys, part_dict in zip([canopy_keys, stem_keys, root_keys], ['canopy_dict', 'stem_dict', 'root_dict']): 
            chap_dict[sp][part_dict] =  {}
            for key in part_keys:
                j = np.where(keys==key)[0][0]+1         # to account for the column of the species
                chap_dict[sp][part_dict][key] = sheet.col_values(spc)[j] 
    return chap_dict

def initialize_plant(sp, params, soil_type):
    # initializing each species with its name (in strings) 
    canopy_dict = params[sp]['canopy_dict']
    stem_dict = params[sp]['stem_dict']
    root_dict = params[sp]['root_dict']
    plant = Whole_plant(species=sp)
    plant.canopy = Canopy(**canopy_dict)
    plant.stem = Stem(**stem_dict)
    plant.soil_root = Soil_root(soil_type=soil_type, **root_dict)
    plant.soil = Soil(soil_type)
    return plant

def get_part(var):
    if var in ['A_canopy', 'Gs_leaf', 'Amax', 'rho', 'lamABA', 'm']: return 'canopy_dict'
    elif var in ['L_stem','A_stem','Ksat_stem','a_stem','c_stem', 'P50_stem']: return 'stem_dict'
    elif var in ['L_root','A_root']: return 'root_dict'
                
def simulate_s_t_norefilling(depths, tRun, dt, sInit, plant, VPD):
    # need to incorporate the case of NO refilling (refilling==False)
    s_t = np.zeros(len(tRun)); s0 = sInit
    assm_t = np.zeros_like(s_t);  assm_t[0] = 0
    px_t = np.zeros_like(s_t)
    pxmin_t = np.zeros_like(s_t)
    E_t = np.zeros_like(s_t)
    gs_t = np.zeros_like(s_t)
    pl_t = np.zeros_like(s_t)
    
    R = plant.canopy.R()
    Ar, Zr, n = plant.soil_root.A_root, plant.soil_root.L_root, plant.soil.n
    Ksat, b = plant.soil.Ksat_soil, plant.soil.b
    
    Px_min = plant.get_fluxes_scalar(VPD,s0,0)[1]
    px_t[0] = Px_min

    for i in range(len(tRun)): 
        R_normed = depths[i]
        Infil_normed = min(R_normed, 1.0-s0) # from the previous step
        
        # Translate s into soil water potential P_soil
        s1 = s0 + Infil_normed  # update with rainfall 
        P_soil = plant.soil.P_soil_solver(s1)
        
        # Calculate fluxes and plant water potentials
        Px = plant.get_fluxes_scalar(VPD,s0,Px_min)[1]
        J = plant.ET_soil_func(P_soil,Px)
        Pl = plant.Pl_ET_func(J, Px, Px_min) 
        
        # calculate the rest of the variables
        ABA = plant.ABA_func(Pl, J)              
        gs = plant.stomata_func(ABA, VPD, Pl) 
             
        # use soil moisture water balance
        L = Ksat*s0**(2*b+3)/(n*Zr*Ar)
        E = J/(n*Zr*Ar)
        s_out = max(s1 - dt*(E + L), 10**(-20)) # drawdown 
        ASM = plant.assimilation_daily(gs)
        
        # check Px to see if it has exceeded Px_min 
        Px_min = min(Px, Px_min)
        
        # update to next step
        s_t[i] = s_out; s0 = s_out
        assm_t[i] = (ASM-R)*dt
        px_t[i] = Px
        pxmin_t[i] = Px_min
        E_t[i] = E
        gs_t[i] = gs
        pl_t[i] = Pl
    # return variables 
    return s_t, assm_t, px_t, pxmin_t, E_t, gs_t, pl_t

def simulate_rainfall(n_trajectories, tRun, dt, lam, gam):
    size = len(tRun)*n_trajectories 
    depthExp = -np.log(1.0-np.random.random(size=size))/gam
    freqUnif = np.random.random(size=size)
    
    depth = np.zeros(size)
    # the occurence of rainfall in any independent interval is lam*dt
    depth[freqUnif<np.tile(lam,size)*dt] = depthExp[freqUnif<np.tile(lam,size)*dt] # rain falls according to prob within an increment
    depth_re = np.reshape(depth, (n_trajectories, len(tRun)))
    return depth_re

def simulate_ps_t_nonlinearized(n_trajectories, tRun, dt, s0, plant, VPD, lam, alpha):
    Ar, Zr, n = plant.soil_root.A_root, plant.soil_root.L_root, plant.soil.n
    gam = (n*Zr*Ar)/(alpha*Ar)

    # set up repositories 
    depth_re = simulate_rainfall(n_trajectories, tRun, dt, lam, gam)
    ps_samples = np.zeros((n_trajectories, len(tRun)))
    assm_samples = np.zeros_like(ps_samples)
    px_samples = np.zeros_like(ps_samples)
    pxmin_samples = np.zeros_like(ps_samples)
    E_samples = np.zeros_like(ps_samples)
    gs_samples = np.zeros_like(ps_samples)
    pl_samples = np.zeros_like(ps_samples)
    for nsim in range(n_trajectories):
        s_t, assm_t, px_t, px_min, E_t, gs_t, pl_t = simulate_s_t_norefilling(depth_re[nsim], tRun, dt, s0, plant, VPD)
        ps_samples[nsim] = s_t
        assm_samples[nsim] = assm_t
        px_samples[nsim] = px_t
        pxmin_samples[nsim] = px_min
        E_samples[nsim] = E_t
        gs_samples[nsim] = gs_t
        pl_samples[nsim] = pl_t
    return ps_samples, assm_samples, px_samples, pxmin_samples, E_samples, gs_samples, pl_samples
