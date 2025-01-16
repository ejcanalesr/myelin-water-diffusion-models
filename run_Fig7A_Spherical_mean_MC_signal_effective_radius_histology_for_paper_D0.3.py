#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Erick Jorge Canales-Rodriguez, 2024

from myelin_water_diffusionMRI_models import *

from   scipy.optimize import minimize_scalar
from   scipy.interpolate import interp1d
from scipy.stats import pearsonr

# Load bvecs
hardi_fname, hardi_bval_fname, hardi_bvec_fname = get_fnames('stanford_hardi')
bvals, bvecs = read_bvals_bvecs(hardi_bval_fname, hardi_bvec_fname)
bvecs = bvecs[bvals>0,:]

bvalues     = np.loadtxt('Data/MC_simulations_Connectome2_waveforms_6bvals/b_val_vector.txt')
BigDelta    = np.loadtxt('Data/MC_simulations_Connectome2_waveforms_6bvals/Delta_vector.txt')
smalldelta  = np.loadtxt('Data/MC_simulations_Connectome2_waveforms_6bvals/delta_vector.txt')
ramptime    = np.loadtxt('Data/MC_simulations_Connectome2_waveforms_6bvals/ramp_vector.txt')

ind_b1  = (bvalues > 0.75) & (bvalues < 0.85)
ind_b2  = (bvalues > 0.95) & (bvalues < 1.05)
ind_b3  = (bvalues > 1.45) & (bvalues < 1.55)
ind_b4 =  (bvalues > 1.95) & (bvalues < 2.05)
ind_b5 =  (bvalues > 2.45) & (bvalues < 2.55)
ind_b6 =  (bvalues > 2.95) & (bvalues < 3.05)

grad_b1 = bvalues[ind_b1]
grad_b2 = bvalues[ind_b2]
grad_b3 = bvalues[ind_b3]
grad_b4 = bvalues[ind_b4]
grad_b5 = bvalues[ind_b5]
grad_b6 = bvalues[ind_b6]

BigDelta1 = BigDelta[ind_b1][0]# ms
BigDelta2 = BigDelta[ind_b2][0]# ms
BigDelta3 = BigDelta[ind_b3][0]# ms
BigDelta4 = BigDelta[ind_b4][0]# ms
BigDelta5 = BigDelta[ind_b5][0]# ms
BigDelta6 = BigDelta[ind_b6][0]# ms

smalldelta1 = smalldelta[ind_b1][0]# ms
smalldelta2 = smalldelta[ind_b2][0]# ms
smalldelta3 = smalldelta[ind_b3][0]# ms
smalldelta4 = smalldelta[ind_b4][0]# ms
smalldelta5 = smalldelta[ind_b5][0]# ms
smalldelta6 = smalldelta[ind_b6][0]# ms

E1 = ramptime[ind_b1][0]# ms
E2 = ramptime[ind_b2][0]# ms
E3 = ramptime[ind_b3][0]# ms
E4 = ramptime[ind_b4][0]# ms
E5 = ramptime[ind_b5][0]# ms
E6 = ramptime[ind_b6][0]# ms

b0      = bvalues == 0

# ---------------------- Define experimental parameters ------------------------
D          = 0.3 # um^2/ms
Nterms     = 40 # For the series
# Lori approximation
Lori        = True
pulse       = 'trapezoid'
#pulse       = 'rectangle'

# -------------------------Generate Signals -----------------------------------#
import scipy.io

#These values are calculated using a fixed grid of radii from 0.1 to 5.0 in steps of 0.1
radius2  = np.linspace(0.1, 5.0, 50)
radius5  = radius2.copy()
radius8  = radius2.copy()
radius10 = radius2.copy()

PM2  = scipy.io.loadmat('Data/Hist/Histvalues_ROI2.mat');  PM2  = PM2 ['Values'][0,:]; PM2  = PM2/np.sum(PM2);
PM5  = scipy.io.loadmat('Data/Hist/Histvalues_ROI5.mat');  PM5  = PM5 ['Values'][0,:]; PM5  = PM5/np.sum(PM5);
PM8  = scipy.io.loadmat('Data/Hist/Histvalues_ROI8.mat');  PM8  = PM8 ['Values'][0,:]; PM8  = PM8/np.sum(PM8);
PM10 = scipy.io.loadmat('Data/Hist/Histvalues_ROI10.mat'); PM10 = PM10['Values'][0,:]; PM10 = PM10/np.sum(PM10);

b          = np.array([0.8, 1.0, 1.5, 2.0, 2.5, 3.0])
Nb         = b.shape[0]
BigDelta   = np.array([BigDelta1, BigDelta2, BigDelta3, BigDelta4, BigDelta5, BigDelta6])
smalldelta = np.array([smalldelta1, smalldelta2, smalldelta3, smalldelta4, smalldelta5, smalldelta6])
E          = np.array([E1, E2, E3, E4, E5, E6])

# ------------------------------------------------------------------------------
# ---------------------------- Load MC signals ---------------------------------
# ------------------------------------------------------------------------------
radius2_discete  = np.round(radius2,1)
radius5_discete  = np.round(radius5,1)
radius8_discete  = np.round(radius8,1)
radius10_discete = np.round(radius10,1)

Signal_sum_MC_2 = 0
MC_2 = np.zeros((Nb, radius2_discete.shape[0]))
for i in range(radius2_discete.shape[0]):
    radius_i  = radius2_discete[i]
    P_i       = PM2[i] # probability

    #---------------------------------------------------------------------------
    Signal_MC_protocol  = np.loadtxt('Data/MC_simulations_Connectome2_waveforms_6bvals/simulationsOuputsD03_G500_WF/myelin_' + str(radius_i) + '0_um_D_0.3_2D_G500_WF_DWI.txt')

    S0_protocol     = np.mean(Signal_MC_protocol[b0])
    SMT_protocol_b1 = np.mean(Signal_MC_protocol[ind_b1])/S0_protocol
    SMT_protocol_b2 = np.mean(Signal_MC_protocol[ind_b2])/S0_protocol
    SMT_protocol_b3 = np.mean(Signal_MC_protocol[ind_b3])/S0_protocol
    SMT_protocol_b4 = np.mean(Signal_MC_protocol[ind_b4])/S0_protocol
    SMT_protocol_b5 = np.mean(Signal_MC_protocol[ind_b5])/S0_protocol
    SMT_protocol_b6 = np.mean(Signal_MC_protocol[ind_b6])/S0_protocol

    Signal_MC = np.array([SMT_protocol_b1, SMT_protocol_b2, SMT_protocol_b3, SMT_protocol_b4, SMT_protocol_b5, SMT_protocol_b6])
    MC_2[:,i] = Signal_MC
    # --------------------------------------------------------------------------

    Signal_sum_MC_2  = Signal_sum_MC_2  + (radius_i*P_i/np.sum(radius2_discete*PM2)) * Signal_MC
#end

Signal_sum_MC_5 = 0
MC_5 = np.zeros((Nb, radius5_discete.shape[0]))
for i in range(radius5_discete.shape[0]):
    radius_i  = radius5_discete[i]
    P_i       = PM5[i] # probability

    #---------------------------------------------------------------------------
    Signal_MC_protocol  = np.loadtxt('Data/MC_simulations_Connectome2_waveforms_6bvals/simulationsOuputsD03_G500_WF/myelin_' + str(radius_i) + '0_um_D_0.3_2D_G500_WF_DWI.txt')

    S0_protocol     = np.mean(Signal_MC_protocol[b0])
    SMT_protocol_b1 = np.mean(Signal_MC_protocol[ind_b1])/S0_protocol
    SMT_protocol_b2 = np.mean(Signal_MC_protocol[ind_b2])/S0_protocol
    SMT_protocol_b3 = np.mean(Signal_MC_protocol[ind_b3])/S0_protocol
    SMT_protocol_b4 = np.mean(Signal_MC_protocol[ind_b4])/S0_protocol
    SMT_protocol_b5 = np.mean(Signal_MC_protocol[ind_b5])/S0_protocol
    SMT_protocol_b6 = np.mean(Signal_MC_protocol[ind_b6])/S0_protocol

    Signal_MC = np.array([SMT_protocol_b1, SMT_protocol_b2, SMT_protocol_b3, SMT_protocol_b4, SMT_protocol_b5, SMT_protocol_b6])
    MC_5[:,i] = Signal_MC
    # --------------------------------------------------------------------------

    Signal_sum_MC_5  = Signal_sum_MC_5  + (radius_i*P_i/np.sum(radius5_discete*PM5)) * Signal_MC
#end

Signal_sum_MC_8 = 0
MC_8 = np.zeros((Nb, radius8_discete.shape[0]))
for i in range(radius8_discete.shape[0]):
    radius_i  = radius8_discete[i]
    P_i       = PM8[i] # probability

    #---------------------------------------------------------------------------
    Signal_MC_protocol  = np.loadtxt('Data/MC_simulations_Connectome2_waveforms_6bvals/simulationsOuputsD03_G500_WF/myelin_' + str(radius_i) + '0_um_D_0.3_2D_G500_WF_DWI.txt')

    S0_protocol     = np.mean(Signal_MC_protocol[b0])
    SMT_protocol_b1 = np.mean(Signal_MC_protocol[ind_b1])/S0_protocol
    SMT_protocol_b2 = np.mean(Signal_MC_protocol[ind_b2])/S0_protocol
    SMT_protocol_b3 = np.mean(Signal_MC_protocol[ind_b3])/S0_protocol
    SMT_protocol_b4 = np.mean(Signal_MC_protocol[ind_b4])/S0_protocol
    SMT_protocol_b5 = np.mean(Signal_MC_protocol[ind_b5])/S0_protocol
    SMT_protocol_b6 = np.mean(Signal_MC_protocol[ind_b6])/S0_protocol

    Signal_MC = np.array([SMT_protocol_b1, SMT_protocol_b2, SMT_protocol_b3, SMT_protocol_b4, SMT_protocol_b5, SMT_protocol_b6])
    MC_8[:,i] = Signal_MC
    # --------------------------------------------------------------------------

    Signal_sum_MC_8  = Signal_sum_MC_8  + (radius_i*P_i/np.sum(radius8_discete*PM8)) * Signal_MC
#end

Signal_sum_MC_10 = 0
MC_10 = np.zeros((Nb, radius10_discete.shape[0]))
for i in range(radius10_discete.shape[0]):
    radius_i  = radius10_discete[i]
    P_i       = PM10[i] # probability

    #---------------------------------------------------------------------------
    Signal_MC_protocol  = np.loadtxt('Data/MC_simulations_Connectome2_waveforms_6bvals/simulationsOuputsD03_G500_WF/myelin_' + str(radius_i) + '0_um_D_0.3_2D_G500_WF_DWI.txt')

    S0_protocol     = np.mean(Signal_MC_protocol[b0])
    SMT_protocol_b1 = np.mean(Signal_MC_protocol[ind_b1])/S0_protocol
    SMT_protocol_b2 = np.mean(Signal_MC_protocol[ind_b2])/S0_protocol
    SMT_protocol_b3 = np.mean(Signal_MC_protocol[ind_b3])/S0_protocol
    SMT_protocol_b4 = np.mean(Signal_MC_protocol[ind_b4])/S0_protocol
    SMT_protocol_b5 = np.mean(Signal_MC_protocol[ind_b5])/S0_protocol
    SMT_protocol_b6 = np.mean(Signal_MC_protocol[ind_b6])/S0_protocol

    Signal_MC = np.array([SMT_protocol_b1, SMT_protocol_b2, SMT_protocol_b3, SMT_protocol_b4, SMT_protocol_b5, SMT_protocol_b6])
    MC_10[:,i] = Signal_MC
    # --------------------------------------------------------------------------

    Signal_sum_MC_10  = Signal_sum_MC_10  + (radius_i*P_i/np.sum(radius10_discete*PM10)) * Signal_MC
#end

#-------------------- Fitting -------------------------------------------------#
from   scipy.optimize import minimize_scalar, fminbound, minimize
from sklearn.linear_model import LinearRegression

def intitial_brute_anat(data_xyz, bvals, D, Nterms, BigDelta, smalldelta, E, Lori, pulse):
    Nb   = bvals.shape[0]
    term = np.zeros((Nb,1))
    Nrad = 500
    x    = np.linspace(0.1, 5.0, Nrad)
    term = np.zeros((Nb, Nrad))
    for i in range(Nb):
        term[i,:] = data_xyz[i] - SMT_signal(Nterms, D, bvals[i], x,  BigDelta[i], smalldelta[i], E[i], Lori, pulse)
    #end i
    cost      = np.sum(term**2, axis=0)
    index_min = np.argmin(cost)
    x_opt     = np.mean(x[index_min])
    return x_opt
#end

def fitting_voxel_anat(data_xyz, bvals, D, Nterms, BigDelta, smalldelta, E, Lori, pulse):
    #x0          = [0.3]
    x0          = intitial_brute_anat(data_xyz, bvals, D, Nterms, BigDelta, smalldelta, E, Lori, pulse)
    bnds        = ( (x0/2, 2*x0),) # bounds
    res         = minimize(obj_fun_anat, x0, method = 'SLSQP', options={'disp': True}, bounds = bnds, args=(data_xyz, bvals, D, Nterms, BigDelta, smalldelta, E, Lori, pulse) )
    reg_sol     = res.x
    return reg_sol[0]
#end fun

def obj_fun_anat(x, data_xyz, bvals, D, Nterms, BigDelta, smalldelta, E, Lori, pulse):
    # Least square error
    Nb   = bvals.shape[0]
    term = np.zeros((Nb,1))
    for i in range(Nb):
        term[i] = data_xyz[i] - SMT_signal(Nterms, D, bvals[i], x,  BigDelta[i], smalldelta[i], E[i], Lori, pulse)
    #end
    cost_fun = np.sum(term**2)
    return cost_fun
#end fun

# ----------------------------------

def intitial_brute_gauss(data_xyz, bvals, D, BigDelta, smalldelta, E, Lori, pulse):
    Nb   = bvals.shape[0]
    term = np.zeros((Nb,1))
    Nrad = 500
    x    = np.linspace(0.1, 5.0, Nrad)
    term = np.zeros((Nb, Nrad))
    for i in range(Nb):
        term[i,:] = data_xyz[i] - SMT_signal_Gaussian(D, bvals[i], x, BigDelta[i], smalldelta[i], E[i], Lori, pulse)
    #end i
    cost      = np.sum(term**2, axis=0)
    index_min = np.argmin(cost)
    x_opt     = np.mean(x[index_min])
    return x_opt
#end

def fitting_voxel_gauss(data_xyz, bvals, D, BigDelta, smalldelta, E, Lori, pulse):
    #               x0
    #x0          = [0.3]
    x0          = intitial_brute_gauss(data_xyz, bvals, D, BigDelta, smalldelta, E, Lori, pulse)
    bnds        = ( (x0/2, 2*x0),) # bounds
    res         = minimize(obj_fun_gauss, x0, method = 'SLSQP', options={'disp': True}, bounds = bnds, args=(data_xyz, bvals, D, BigDelta, smalldelta, E, Lori, pulse)  )
    reg_sol     = res.x
    return reg_sol[0]
#end fun

def obj_fun_gauss(x, data_xyz, bvals, D, BigDelta, smalldelta, E, Lori, pulse):
    # Least square error
    Nb   = bvals.shape[0]
    term = np.zeros((Nb,1))
    for i in range(Nb):
        term[i] = data_xyz[i] - SMT_signal_Gaussian(D, bvals[i], x, BigDelta[i], smalldelta[i], E[i], Lori, pulse)
    #end
    cost_fun = np.sum(term**2)
    return cost_fun
#end fun

# ----------------------------------

def fitting_voxel_MC(Signal_sum_MC, MC, radius_discrete):
    MSE = np.sum((MC - np.array(Signal_sum_MC, ndmin=2).T)**2, axis=0)
    # --- Using brute force, disretization on the grid
    # radius = radius_discrete[np.argmin(MSE)]

    # --- Using Cubic B-spline interpolation
    f2  = interp1d(radius_discrete, MSE, kind='cubic')
    res = minimize_scalar(f2, method='Bounded', bounds=(0.2, 4.9))
    radius = res.x
    return radius
#end

print('Fitting the MC signal to a single cylinder using the MC signals')
Est_effective_MC_radius2  = fitting_voxel_MC(Signal_sum_MC_2,  MC_2,  radius2_discete)
Est_effective_MC_radius5  = fitting_voxel_MC(Signal_sum_MC_5,  MC_5,  radius5_discete)
Est_effective_MC_radius8  = fitting_voxel_MC(Signal_sum_MC_8,  MC_8,  radius8_discete)
Est_effective_MC_radius10 = fitting_voxel_MC(Signal_sum_MC_10, MC_10, radius10_discete)
effective_radius_MC       = np.array([Est_effective_MC_radius2, Est_effective_MC_radius5, Est_effective_MC_radius8, Est_effective_MC_radius10])

print('Fitting the MC signal to a single cylinder using the General model')
# General model
Est_effective_anat_radius2  = fitting_voxel_anat(Signal_sum_MC_2,  b, D, Nterms, BigDelta, smalldelta, E, Lori, pulse)
Est_effective_anat_radius5  = fitting_voxel_anat(Signal_sum_MC_5,  b, D, Nterms, BigDelta, smalldelta, E, Lori, pulse)
Est_effective_anat_radius8  = fitting_voxel_anat(Signal_sum_MC_8,  b, D, Nterms, BigDelta, smalldelta, E, Lori, pulse)
Est_effective_anat_radius10 = fitting_voxel_anat(Signal_sum_MC_10, b, D, Nterms, BigDelta, smalldelta, E, Lori, pulse)
effective_radius_anat       = np.array([Est_effective_anat_radius2, Est_effective_anat_radius5, Est_effective_anat_radius8, Est_effective_anat_radius10])

print('Fitting the MC signal to a single cylinder using the Gaussian approximation')
# Gaussian model
Est_effective_gauss_radius2  = fitting_voxel_gauss(Signal_sum_MC_2,  b, D, BigDelta, smalldelta, E, Lori, pulse)
Est_effective_gauss_radius5  = fitting_voxel_gauss(Signal_sum_MC_5,  b, D, BigDelta, smalldelta, E, Lori, pulse)
Est_effective_gauss_radius8  = fitting_voxel_gauss(Signal_sum_MC_8,  b, D, BigDelta, smalldelta, E, Lori, pulse)
Est_effective_gauss_radius10 = fitting_voxel_gauss(Signal_sum_MC_10, b, D, BigDelta, smalldelta, E, Lori, pulse)
effective_radius_gauss       = np.array([Est_effective_gauss_radius2, Est_effective_gauss_radius5, Est_effective_gauss_radius8, Est_effective_gauss_radius10])

# Histology
mean_radius   = np.array([np.sum(radius2_discete*PM2), np.sum(radius5_discete*PM5), np.sum(radius8_discete*PM8), np.sum(radius10_discete*PM10)])
mean_radius21 = np.array([np.sum((radius2_discete**2)*PM2), np.sum((radius5_discete**2)*PM5), np.sum((radius8_discete**2)*PM8), np.sum((radius10_discete**2)*PM10)])/mean_radius
mean_radius31 = np.sqrt(np.array([np.sum((radius2_discete**3)*PM2), np.sum((radius5_discete**3)*PM5), np.sum((radius8_discete**3)*PM8), np.sum((radius10_discete**3)*PM10)])/mean_radius)

# Write file
import csv
with open('effective_radius_MC_signals_D_' + str(D) + '.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter = ',', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerow(['effective_radius_MC', effective_radius_MC])
    spamwriter.writerow(['effective_radius_anat', effective_radius_anat])
    spamwriter.writerow(['effective_radius_gauss', effective_radius_gauss])
    spamwriter.writerow(['mean_radius', mean_radius])
    spamwriter.writerow(['mean_radius21', mean_radius21])
    spamwriter.writerow(['mean_radius31', mean_radius31])

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
#
ax1 = plt.figure(figsize=(4.8,4.8))

res = pearsonr(mean_radius, effective_radius_anat)
plt.plot(mean_radius, effective_radius_anat, 'o', label='$a_{eff}$' + ' vs ' + r"$<a>$" + r",    $PCC$ = " + str(np.round(res.statistic,3)) + r" ($p$ = " + str(np.round(res.pvalue,3)) + r")", color=colors[0])

X, Y = mean_radius.reshape(-1,1), effective_radius_anat.reshape(-1,1)
plt.plot( X, LinearRegression().fit(X, Y).predict(X) , color=colors[0])

res = pearsonr(mean_radius21, effective_radius_anat)
plt.plot(mean_radius21, effective_radius_anat, 'o', label='$a_{eff}$' + ' vs ' + r"$\frac{<a^2>}{<a>}$" + r",    $PCC$ = " + str(np.round(res.statistic,3)) + r" ($p$ = " + str(np.round(res.pvalue,3)) + r")", color=colors[1])

X1, Y1 = mean_radius21.reshape(-1,1), effective_radius_anat.reshape(-1,1)
plt.plot( X1, LinearRegression().fit(X1, Y1).predict(X1), color=colors[1] )

res = pearsonr(mean_radius31, effective_radius_anat)
plt.plot(mean_radius31, effective_radius_anat, 'o', label='$a_{eff}$' + ' vs ' + r"$(\frac{<a^3>}{<a>})^{\frac{1}{2}}$" + r",    $PCC$ = " + str(np.round(res.statistic,3)) + r" ($p$ = " + str(np.round(res.pvalue,3)) + r")", color=colors[2])
X2, Y2 = mean_radius31.reshape(-1,1), effective_radius_anat.reshape(-1,1)
plt.plot( X2, LinearRegression().fit(X2, Y2).predict(X2), color=colors[2] )

plt.xlabel(r"$Moments$" + "-" + "$based$" + " " + "$radius$" + " $(\mu m)$", fontsize=14)
plt.ylabel(r"$a_{eff}$" + " $(\mu m)$", fontsize=14)

ax = plt.gca()
ax.legend(fontsize=12, ncol=1, handleheight=0, labelspacing=0.75, frameon=False)
ax.set_ylim([0.0, 1.2])
ax.set_xlim([0.0, 1.2])

plt.title(r"Effective vs Moments-based radii " + r"($D$ = " + str(D) + "$\mu m^2/ms$)", fontsize=14)
plt.savefig('Fig7A_Effective_vs_Moments_based_radii_using_MC_signal_D_'  + str(D)  + '.png', bbox_inches='tight', dpi=300)

#ax2 = plt.figure(figsize=(6,6))
ax2 = plt.figure(figsize=(4.8,4.8))

plt.plot(mean_radius, effective_radius_MC, 'o', label='$a_{eff-MC}$' + ' vs ' + r"$<a>$" + r",    $PCC$ = " + str(np.round(np.corrcoef(mean_radius, effective_radius_MC)[0,1],3)), color=colors[0])
X, Y = mean_radius.reshape(-1,1), effective_radius_MC.reshape(-1,1)
plt.plot( X, LinearRegression().fit(X, Y).predict(X) , color=colors[0])

plt.plot(mean_radius21, effective_radius_MC, 'o', label='$a_{eff-MC}$' + ' vs ' + r"$\frac{<a^2>}{<ra>}$" + r",   $PCC$ = " + str(np.round(np.corrcoef(mean_radius21, effective_radius_MC)[0,1],3)), color=colors[1])
X1, Y1 = mean_radius21.reshape(-1,1), effective_radius_MC.reshape(-1,1)
plt.plot( X1, LinearRegression().fit(X1, Y1).predict(X1), color=colors[1] )

plt.plot(mean_radius31, effective_radius_MC, 'o', label='$a_{eff-MC}$' + ' vs ' + r"$(\frac{<a^3>}{<a>})^{\frac{1}{2}}$" + r", $PCC$ = " + str(np.round(np.corrcoef(mean_radius31, effective_radius_MC)[0,1],3)), color=colors[2])
X2, Y2 = mean_radius31.reshape(-1,1), effective_radius_MC.reshape(-1,1)
plt.plot( X2, LinearRegression().fit(X2, Y2).predict(X2), color=colors[2] )

plt.xlabel(r"$radius$" + " $(\mu m)$", fontsize=14)
plt.ylabel(r"$a_{eff-MC}$" + " $(\mu m)$", fontsize=14)

ax = plt.gca()
ax.legend(fontsize=12, ncol=1, handleheight=0, labelspacing=0.75, frameon=False)
ax.set_ylim([0.0, 1.2])
ax.set_xlim([0.0, 1.2])

plt.title(r"Effective vs Moments-based radii " + r"($D$ = " + str(D) + "$\mu m^2/ms$)", fontsize=14)
#plt.savefig('Supp_Fig_7_Effective_vs_Moments_based_radii_using_MC_signal_and_MC_dictionary_D_'  + str(D) + '.png', bbox_inches='tight', dpi=300)

ax3 = plt.figure(figsize=(4.8,4.8))
plt.plot(effective_radius_anat, effective_radius_gauss, 'o', label=r"$<radius>$ " + " General vs Gaussian aprox", color=colors[3])
X, Y = effective_radius_anat.reshape(-1,1), effective_radius_gauss.reshape(-1,1)
plt.ylabel(r"$a_{eff} (Gaussian approximation)$" + " $(\mu m)$", fontsize=14)
plt.xlabel(r"$a_{eff} (General equation)$" + " $(\mu m)$", fontsize=14)
plt.plot( X, LinearRegression().fit(X, Y).predict(X) , color=colors[3])

plt.show()
