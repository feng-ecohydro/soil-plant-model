'''

This model couples stomatal conductance, plant hydraulics, ABA production, and soil water balance to quantify 
(1) the hydraulic risk and (2) net cumulative carbon assimilation of a plant over a designated dry period.
It does not allow for xylem refilling. For more details of model implementation, please refer to: 

Feng et al. (2018), The ecohydrological context of drought and classification of plant responses, Ecology Letters.

This code was written in Python 2.7 and requires the Python packages numpy, scipy, matplotlib, and xlrd.
They can be downloaded through an Anaconda Distribution, which also includes the development environment Spyder. 
See: https://www.anaconda.com/download/ and https://docs.anaconda.com/anaconda/packages/pkg-docs.html (for a list of packages). 

Author contact: feng@umn.edu
Date: 2018/07/14

'''
import numpy as np
import scipy.optimize as opt
import params_constants 
from params_soil import soil_dict

# Prerequisite thresholds, parameters, functions
low_val_pos = 10**(-20)     # low positive values used to constrain stomatal conductance gs
pxf = 2.0                   # a multiplier to P50 beyond which xylem water potential is not allowed 
pl_min  = -50.0             # minimum allowable leaf water potential
# function for clipping an input value between designated minimum and maximum
clamp = lambda minimum, x, maximum: max(minimum, min(x, maximum))

class Canopy: 
    def __init__(self, Gs_leaf, A_canopy, Amax, rho, lamABA, m): 
        self.Gs_leaf = Gs_leaf                      # max leaf conductance for H2O, mol/m2/s
        self.A_canopy = A_canopy                    # canopy area, m2
        self.Amax = Amax                            # maximum daily assimilation, mol/m2/d
        self.rho = rho                              # ratio of Rmax/Amax
        self.m = m                                  # exponent on the stomatal response curve
        self.lamABA = lamABA                        # ABA production rates, mol ABA/m2-s-MPa
        
    def g_canopy(self, gs):
        # Maximum canopy conductance per day, m3/d
        Sd, rhoH2O, Mw = params_constants.Sd, params_constants.rhoH2O, params_constants.Mw
        return gs*self.A_canopy*Mw*Sd/(rhoH2O)
    
    def R(self):
        # Daily respiration rate, as constant ratio rho of Amax, m3/d
        return self.rho*self.Amax
    
    def carboxyl_constant(self):
        # Carboxylation capacity, k; mol/m2-s
        # Amax should be inputed as mol/m2/d
        # CO2 concentration = 400 ppm = 400E-6 mol/mol
        Sd, CO2conc = params_constants.Sd, params_constants.CO2conc
        gmax, Amax, R = self.Gs_leaf, self.Amax/Sd, self.R()/Sd # Amax and R converted to per second
        return -Amax*gmax/(1.6 * (Amax-R) - gmax*CO2conc)
    
    def kmax_canopy(self):
        # Carboxylation capacity converted to daily canopy level; mol/d
        kmax = self.carboxyl_constant()
        Sd, rhoH2O, Mw = params_constants.Sd, params_constants.rhoH2O, params_constants.Mw
        return kmax*self.A_canopy*Mw*Sd/(rhoH2O)
    
class Stem: 
    def __init__(self, L_stem, A_stem, Ksat_stem, P50_stem, a_stem=None, c_stem=None, plc_form=None, stem_safety=None):
        self.L_stem = L_stem            # length of conducting stem, m 
        self.A_stem = A_stem            # sapwood area, m2
        self.Ksat_stem = Ksat_stem      # maximum hydraulic conductivity, kg/m2-s-MPa 
        self.P50_stem = P50_stem        # water potential at 50% loss in maximum xylem conductivity, MPa
        self.plc_form = plc_form        # functional form of xylem vulnerability, either 'exp' for exponential or 'plc' for exponential sigmoidal
        self.a_stem = a_stem            # parameter of xylem vulnerability in 'plc' form
        self.c_stem = c_stem            # parameter of xylem vulnerability in 'exp' form
    
    def kmax_stem(self):
        # maximum hydraulic conductivity in the stem, m3/d
        Sd, rhoH2O = params_constants.Sd, params_constants.rhoH2O
        return self.Ksat_stem*Sd*self.A_stem/(self.L_stem*rhoH2O)
    
    def k_stem(self, P_stem):
        # Stem conductivity, m3/d-MPa
        kmax, a, c, P50 = self.kmax_stem(), self.a_stem, self.c_stem, self.P50_stem
        if c:  # based on exponential function 
            return kmax*np.exp( -c * clamp(pxf*P50, P_stem,0) )
        else: 
            # based on exponential-sigmoidal function 
            # constrain the domain of P to be negative
            return kmax*(1.0-1.0/(1.0+np.exp( a*( clamp(pxf*P50, P_stem,0) - P50) )))
        
    def k_stem_Kirchhoff(self, P):
        # Kirchoff transform - integrate vulnerability curve from -infinity to P_stem (P) 
        kmax, a, c, P50 = self.kmax_stem(), self.a_stem, self.c_stem, self.P50_stem
        if c: # based on exponential function 
            k_stem = kmax*( np.exp(-c* clamp(pxf*P50,P,0) )/-c )
        else: 
            # based on exponential-sigmoidal function 
            # constrain the domain of P to be negative
            k_stem = kmax*(np.log( np.exp(-a*P50) + np.exp(- a* clamp(pxf*P50,P,0)) )/a + P)
        return k_stem
    
    def Px_solvent(self, P_stem, percent):
        # helper function for solving for P_stem corresponding to designated percent loss of conductivity
        return self.k_stem(P_stem) - (1.0-percent)*self.kmax_stem()
        
class Soil_root: 
    def __init__(self, L_root, A_root, d_root, soil_type='loamy_sand'):
        # select soil types: "sand", "loamy_sand", "sandy_loam", "loam", "clay" to set hydrological properties 
        self.soil_type = soil_type 
        self.L_root = L_root        # mean rooting depth, Zr, m
        self.A_root = A_root        # rooting area, m2
        self.d_root = d_root        # fine root diameter, m
        
    def kmax_root(self):
        # Maximum soil-root hydraulic conductivity, m3/d-MPa
        Ksat_soil = soil_dict[self.soil_type]['Ksat']
        nu, rhoH2O, g = params_constants.nu, params_constants.rhoH2O, params_constants.g
        return nu*Ksat_soil/(rhoH2O*g)*np.sqrt(self.A_root/(self.d_root*self.L_root))
    
    def k_soil(self, P_soil):
        # Soil-root conductivity based on soil water potential, m3/d-MPa
        b = soil_dict[self.soil_type]['b']
        Psat_soil = soil_dict[self.soil_type]['Ps']
        kmax = self.kmax_root()
        c = 2.0*b+3.0; d =4.0
        # Cap plant water potential by soil water potential
        Pcorrected = min(Psat_soil, P_soil) 
        return kmax*(Psat_soil/Pcorrected)**((c-d)/b)
    
    def k_soil_Kirchhoff(self, P):
        # Kirchhoff transform - integrated k_soil from -inf to P (which cannot be greater than Psat_soil)
        b = soil_dict[self.soil_type]['b']
        Psat_soil = soil_dict[self.soil_type]['Ps']
        kmax = self.kmax_root()
        c = 2.0*b+3.0; d =4.0
        # Cap plant water potential by soil water potential
        Pcorrected = min(Psat_soil, P) 
        return kmax*(b/(b-c+d))*Pcorrected*(Psat_soil/Pcorrected)**((c-d)/b)
    
class Soil:
    def __init__(self, soil_type):
        # select soil types: "sand", "loamy_sand", "sandy_loam", "loam", "clay" to set hydrological properties 
        self.soil_type = soil_type 
        self.Ksat_soil = soil_dict[self.soil_type]['Ksat']  # Saturated soil hydraulic conductivity, m/d
        self.Psat_soil = soil_dict[self.soil_type]['Ps']    # Soil wter potential near soil saturation, MPa
        self.b = soil_dict[self.soil_type]['b']             # Exponent on soil water retention curve
        self.n = soil_dict[self.soil_type]['n']             # Soil porosity
        self.sh = soil_dict[self.soil_type]['sh']           # Hygroscopic point
        self.sfc = soil_dict[self.soil_type]['sfc']         # Field capacity
        
    def s_solver(self, P_soil):
        # returns soil moisture based on soil water potential 
        return (P_soil/self.Psat_soil)**(-1.0/self.b)
    
    def P_soil_solver(self, s):
        # return soil water potential based on soil moisture
        return self.Psat_soil*s**(-self.b)
    
class Whole_plant: 
    def __init__(self, species):
        self.species = species       # Species label
        self.a_seq = 10**(-4)        # ABA sequestration rate, molH2O/m2-s
        self.D0 = 0.0167             # Parameter combination, dimensionless, see Dewar (2002)
        self.beta = 1.48*10**(-4)    # ABA sensitivity parameter, m3/mol ABA
        self.Vw = 18*10**(-6)        # Partial molal volume of water m3/mol
    
    @property
    def canopy(self):
        return self._canopy
    @canopy.setter
    def canopy(self, canopy_object):
        self._canopy = canopy_object
        
    @property
    def stem(self):
        return self._stem
    @stem.setter
    def stem(self, stem_object):
        self._stem =  stem_object
        
    @property
    def soil_root(self):
        return self._soil_root
    @soil_root.setter
    def soil_root(self, soil_root_object):
        self._soil_root =  soil_root_object
        
    @property
    def soil(self):
        return self._soil
    @soil.setter
    def soil(self, soil_object):
        self._soil =  soil_object
    
    def ABA_func(self, P_leaf, ET):
        # Leaf xylem ABA concentration - Equation (6) in Feng (2018)
        lamABA = self.canopy.lamABA     # ABA production rate
        a_seq = self.a_seq              # ABA sequestration rate
        Vw = self.Vw                    # Molal volume of water
        Sd = params_constants.Sd        # Seconds in day
        A_canopy = self.canopy.A_canopy # Canopy area
        return  -lamABA*P_leaf/ ( ET/(A_canopy*Sd) + a_seq*Vw ) # molABA/m2-s
    
    def stomata_func(self, ABA, VPD, P_leaf):
        # Stomatal conductance of CO2 - Equation (5) in Feng (2018)
        rho, kcarbx, gmax = self.canopy.rho, self.canopy.carboxyl_constant(), self.canopy.Gs_leaf
        # a1 is a constant used to designate gmax based on Amax, through effects on the carboxylation capacity
        # gs should be maximized at gmax if a1 is calculated by setting Ds=D0
        a1 = 2.0*gmax/(1.6*kcarbx*(1.0+rho)) 
        D0, beta = self.D0, self.beta
        P_baro = params_constants.P_baro
        # variables are restricted within reasonable ranges
        return a1*(1.0+rho)*kcarbx/(1.0+(VPD/(P_baro*D0))) * np.exp( - beta*max(ABA,0))* np.exp(self.canopy.m * max(P_leaf,pl_min))
    
    def log_stomata_func(self, ABA, VPD, P_leaf):
        # Logarithm of the previous stomatal conductance function (Equation (5) in Feng (2018))
        rho, kcarbx, gmax = self.canopy.rho, self.canopy.carboxyl_constant(), self.canopy.Gs_leaf
        a1 = 2.0*gmax/(1.6*kcarbx*(1.0+rho)) 
        D0, beta = self.D0, self.beta
        P_baro = params_constants.P_baro
        return np.log( a1*(1.0+rho)*kcarbx/(1.0+(VPD/(P_baro*D0))) ) - max(ABA,0)*beta + self.canopy.m * max(P_leaf,pl_min) 
    
    def ET_soil_func(self, P_soil, P_stem):
        # Transpiration from soil to stem, m3/d; Equation (1) of Feng (2018)
        k_soilKirch, P50 = self.soil_root.k_soil_Kirchhoff, self.stem.P50_stem 
        # Again, constrain variables within reasonable range of values 
        return k_soilKirch( min(P_soil,0) ) - k_soilKirch( clamp(pxf*P50,P_stem,0) ) 
    
    def ET_stem_func(self, P_stem, P_leaf, Px_min):
        # Transpiration from stem to canopy, m3/d; Equation (1) of Feng (2018)
        k_stemKirch = self.stem.k_stem_Kirchhoff 
        k_stem, P50 = self.stem.k_stem, self.stem.P50_stem
        Px_clamped = clamp(pxf*P50, P_stem, 0)
        # To account for refilling
        return k_stemKirch(Px_clamped) - k_stemKirch(min(P_leaf, P_stem)) + max(0, k_stem(Px_min)*(P_stem - Px_min))
    
    def ET_canopy_func(self, gs, VPD):
        # Transpiration from the canopy; m3/d
        g_canopy = self.canopy.g_canopy
        P_baro = params_constants.P_baro
        return  g_canopy(gs)*(VPD/P_baro) 
    
    def Pl_ET_func(self, ET, P_stem, Px_min):
        # Solving for leaf water potential Pl given (ET, P_stem, and Px_min)
        # Px_min is the minimum stem water potential experienced by the plant 
        a, c, P50, kmax, k_stem = self.stem.a_stem, self.stem.c_stem, self.stem.P50_stem, self.stem.kmax_stem(), self.stem.k_stem
        
        # change the effective ET and P_stem processed over the untraversed area of the vulnerability curve (v.c.) to account for no refilling 
        ET = ET - max(0, k_stem(Px_min)*(P_stem - Px_min))  # effective ET; if P_stem > Px_min, then account for area over kmax
        P_stem = min(P_stem, Px_min)                        # effective P_stem; if P_stem > Px_min, then start traversing v.c. from Px_min
        
        if c: # calculate leaf water potential from vulnerability curve ("exp" form)
            Pl = 1.0/c * np.log(1.0/(ET/kmax*c + np.exp(-c*P_stem)))
        else: # vulnerability curve in "plc" form
            h = 1.0 - np.exp(a*ET/kmax) + np.exp(a*(P_stem - P50))
            if h <=0: Pl = 1.1*P_stem
            else: Pl = -ET/kmax + P50 + 1.0/a * np.log(h)   
        return Pl
        
    def assimilation_daily(self,gs):
        # Returns daily assimilation rate based on carboxylation capacity, mol/m2-d
        Sd, CO2conc, kcarbx = params_constants.Sd, params_constants.CO2conc, self.canopy.carboxyl_constant()
        R = self.canopy.R()/Sd
        A = kcarbx *(gs*CO2conc + 1.6*R) /(1.6*kcarbx + gs) * Sd
        return A
    
    def flux_solvent(self, (ET, P_stem, P_leaf, gs, ABA), P_soil, VPD, Px_min):
        # Helper function for solving coupled soil-plant hydraulic system 
        # A can be solved later as function of gs - Equation (9) of Feng (2018)
        # A - gs*kcarbx *CO2conc/(1.6*kcarbx + gs)     # mol/m2-s 
        # Constraints: gs need to be greatern than 0, minimum P_stem
        
        return (ABA - self.ABA_func(P_leaf, ET), 
                np.log( max(gs, low_val_pos) ) - self.log_stomata_func(ABA, VPD, P_leaf),
                ET - self.ET_soil_func(P_soil, P_stem), 
                ET - self.ET_stem_func(P_stem, P_leaf, Px_min), 
                ET - self.ET_canopy_func(gs, VPD)  )

    def flux_solver(self, P_soil, VPD, Px_min, init):
        # Solves coupled soil-plant hydraulic system - system of 3 equations describing internal water states (ET, P_stem, P_leaf) 
        # based on an input of soil water potential and VPD
        sol = opt.root(self.flux_solvent, init, (P_soil, VPD, Px_min))
        return sol
    
    def flux_solver_s(self, s, VPD, Px_min, init):
        # Solves coupled soil-plant hydraulic system - system of 3 equations describing internal water states (ET, P_stem, P_leaf) 
        # based on an input of soil moisture and VPD
        P_soil = self.soil.P_soil_solver(s)
        sol =  self.flux_solver(P_soil, VPD, Px_min, init)
        return sol
    
    def calculate_from_Px(self, Px, VPD, Px_min, Ps):
        # Calculates system variables based on designated stem water potential Px. 
        J = self.ET_soil_func(Ps,Px)            # Transpiration from soil to stem
        Pl = self.Pl_ET_func(J, Px, Px_min)     # Leaf water potential
        ABA = self.ABA_func(Pl, J)              # ABA concentration
        gs = self.stomata_func(ABA, VPD, Pl)    # Stomatal conductance  
        ET = self.ET_canopy_func(gs, VPD)       # Transpiration from canopy
        return J, Pl, ABA, gs, ET
    
    def get_init_conditions(self, VPD, si, Px_min, Px_arr):
        # Outputs guesses for initial conditions used to solve the system, based on inputs of soil moisture, VPD, and minimum stem water potential Px_min
        # Px_arr is the range of stem water potentials over which potential initial conditions are tested. 
        Ps = self.soil.P_soil_solver(si)
        Inits_arr = np.zeros((len(Px_arr),5)); Inits_arr2 = np.zeros_like(Inits_arr)
        diffET_arr = np.zeros_like(Px_arr); J_arr = np.zeros_like(Px_arr); ET_arr = np.zeros_like(Px_arr)
        for i, Px in enumerate(Px_arr): 
            J, Pl, ABA, gs, ET = self.calculate_from_Px(Px, VPD, Px_min, Ps)
            Inits_arr[i] = [ET, Px, Pl, gs, ABA]
            Inits_arr2[i] = [J, Px, Pl, gs, ABA]
            J_arr[i] = J
            ET_arr[i] = ET
            diffET_arr[i] = J-ET
        # find when difference between water supply and demand is minimal
        j = np.argmin(np.abs(diffET_arr))
        init_list = []
        try: 
            for k in [0,1,-1]:
                init_list.extend([Inits_arr[j+k]])
                init_list.extend([Inits_arr2[j+k]])
        except IndexError:  # if minimum difference for water supply and demand is found close to the edge of Px range
            init_list.extend([Inits_arr[-1], Inits_arr2[-1]]) 
            init_list.extend([Inits_arr[0], Inits_arr2[0]])
        return init_list
    
    def get_fluxes(self, VPD, s, Px_min, inits=None):
        # returns solutions to the soil-plant system (E, P_stem, P_leaf) along with correponsding P_soil
        Flux = np.zeros((len(s), 6))
        P_soil = np.zeros(len(s))
        for i, si in enumerate(s):
            P_soil[i] = self.soil.P_soil_solver(si) 
            Flux[i] = self.get_fluxes_scalar(VPD, si, Px_min, inits)
        return Flux, P_soil
    
    def get_fluxes_scalar(self, VPD, s, Px_min, Flux_inits=None):
        # returns solutions to the soil-plant system (E, P_stem, P_leaf) along with correponsding P_soil
        
        # Set a range of stem water potentials Px over which initial conditions can be guessed
        Px_arr = self.soil.Psat_soil*np.linspace(s, self.soil.sh-0.001, 50)**(-self.soil.b)
        
        sol_flag = False # flag will be turned to True only when solution is found with inputted guesses
        if Flux_inits is not None: # if initial conditions have already been included in the arguments
            sol0 = self.flux_solver_s(s, VPD, Px_min, Flux_inits)
            sol_flag = sol0.success
            if sol_flag: 
                ET, P_stem, P_leaf, gs, ABA = sol0.x
        # if initial conditions have not been designated, or if previous attempt at solving was not successful, 
        # then we used guessed initial conditions to solve the system
        if (Flux_inits is None) or (sol_flag is False):  
            Flux_inits_generated = self.get_init_conditions(VPD, s, Px_min, Px_arr)
            sol1 = self.flux_solver_s(s, VPD, Px_min, Flux_inits_generated[0])

            if sol1.success: 
                ET, P_stem, P_leaf, gs, ABA = sol1.x
            else: 
                for _, init in enumerate(Flux_inits_generated[1:]):
                    sol = self.flux_solver_s(s, VPD, Px_min, init) 
                    if sol.success: 
                        ET, P_stem, P_leaf, gs, ABA = sol.x
                        break 
                    # the default is the best initial guess
                    ET, P_stem, P_leaf, gs, ABA = Flux_inits_generated[0]
        # calculate assimilation independently
        A = self.assimilation_daily(gs)
        # return variables 
        return [ET, P_stem, P_leaf, gs, ABA, A]
