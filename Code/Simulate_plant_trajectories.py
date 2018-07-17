from utility_functions import import_traits_data, initialize_plant, simulate_ps_t_nonlinearized
from params_soil import soil_dict
import matplotlib.pyplot as plt
import numpy as np

# set soil conditions 
soil_type = 'loamy_sand'
smin = soil_dict[soil_type]['sh']
sfc = soil_dict[soil_type]['sfc']
sst = soil_dict[soil_type]['sst']
n = soil_dict[soil_type]['n']
Ks = soil_dict[soil_type]['Ksat']

# set species traits
traits = import_traits_data()
juni_plant = initialize_plant('JUNI', traits, soil_type)
pine_plant = initialize_plant('PINE', traits, soil_type)

def plot_trajectories(plant, tmax, VPD, newfig=False, lam=0.1, alpha=0.01, s0=sst, dt=0.1, n_trajectories=100): 
    if newfig: plt.figure(figsize=(6,7))
    tRun =  np.arange(0,tmax+dt,dt)
    Amax = plant.canopy.Amax
    R = plant.canopy.R()
    
    # simulate trajectories
    ps, assm, _, pxmin, _, _, _ = simulate_ps_t_nonlinearized(n_trajectories, tRun, dt, s0, plant, VPD, lam, alpha)
    ps_mean = np.mean(ps, axis=0)
    assm_cumsum = np.cumsum(assm, axis=1)/((Amax-R)*np.shape(assm)[1]*dt)
    pxmin_mean = np.mean(pxmin, axis=0)
 
    plt.figure(figsize=(5,8))
    plt.subplot(3,1,1)
    for i in range(len(ps)):
        plt.plot(tRun, ps[i], lw=0.4, color='darkgray')
    plt.plot(tRun, ps[0], lw=0.5, color='blue')
    plt.plot(tRun, ps_mean, lw=1.5, color='blue') 
    plt.ylabel('Mean soil moisture')
    
    plt.subplot(3,1,2)
    for i in range(len(ps)):
        plt.plot(tRun, pxmin[i], lw=0.4, color='darkgray')
    plt.plot(tRun, pxmin[0], lw=0.5, color='red')
    plt.plot(tRun, pxmin_mean, lw=1.5, color='red')
    plt.ylabel('Minimum stem \n water potential (MPa)')
    plt.ylim(-3,0)
    
    plt.subplot(3,1,3)
    for i in range(len(ps)):
        plt.plot(tRun, assm_cumsum[i], lw=0.4, color='darkgray')
    plt.plot(tRun, np.mean(assm_cumsum,axis=0), lw=1.5, color='black')
    plt.ylabel('Net cumulative \n assimilation')
    plt.xlabel('Days')
    
    # output HR and CA
    print plant.species
    pxCrit = plant.stem.P50_stem
    print 'HR: ', len(pxmin[pxmin<pxCrit])/float(np.shape(pxmin)[0]*np.shape(pxmin)[1])
    print 'CA: ', np.mean(assm_cumsum,axis=0)[-1]
    plt.tight_layout()

if __name__ == '__main__':
    plot_trajectories(pine_plant, tmax=90, VPD=2.0, n_trajectories=20)
    plot_trajectories(juni_plant, tmax=90, VPD=2.0, n_trajectories=20)
    plt.show()
    
