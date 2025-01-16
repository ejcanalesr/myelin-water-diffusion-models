#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Erick Jorge Canales-Rodriguez, 2024

from myelin_water_diffusionMRI_models import *

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

# ---------------------- Define experimental parameters ------------------------
#D          = 0.8
D          = 0.3
Nterms     = 40 # For the series
# Lori approximation
Lori        = True
pulse       = 'trapezoid'
#pulse       = 'rectangle'

# -------------------------Generate Signals -----------------------------------#
radius   = np.array([0.3, 1.0, 3.0])
Nb       = 201
b        = np.linspace(0, 100, Nb)

Signal_anat = np.zeros((Nb, radius.shape[0]))
Rad_signal  = np.zeros_like(Signal_anat)
for i in range(Nb):
    print(i)
    bi = b[i]
    Signal_anat[i,:]   = SMT_signal(Nterms, D, b[i], radius, BigDelta6, smalldelta6, E6, Lori, pulse)
#end

# ------------------------ Plot results ----------------------------------------
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
import matplotlib.patches as patches
#
ax1 = plt.figure(figsize=(7.5,5))

plt.plot(b, Signal_anat[:,0], label=r"$radius$ = " + str(round(radius[0],1)) +  " $\mu m$", color=colors[0])
plt.plot(b, Signal_anat[:,1], label=r"$radius$ = " + str(round(radius[1],1)) +  " $\mu m$", color=colors[1])
plt.plot(b, Signal_anat[:,2], label=r"$radius$ = " + str(round(radius[2],1)) +  " $\mu m$", color=colors[2])

ax = plt.gca()
ax.legend(fontsize=12, ncol=1, handleheight=0, labelspacing=0.75, frameon=False, title=r"$\Delta$ = " + str(round(BigDelta6,3)) + " ms, $\delta$ = " + str(round(smalldelta6,3)) + " ms, " + "$D$ = " + str(D) + "$\mu m^2/ms$", title_fontsize=12)

ax.set_yscale('log')

plt.xticks(fontsize=11.5)
plt.yticks(fontsize=11.5)
plt.ylabel(r"Log of Normalized Signal Amplitude", fontsize=13)
plt.xlabel(r"$b$" + " $(ms/ \mu m^2)$", fontsize=14)
#
bmax = 3.0

plt.title(r"Myelin Water Spherical Mean dMRI Signal", fontsize=14)
rect = patches.Rectangle((-0.2, 0.1), bmax*1.1, 0.95, linewidth=1, edgecolor='k', facecolor='none')
ax.add_patch(rect)

ax2 = plt.axes([.6565, .635, .24, .24])
ind = b<=bmax
ax2.plot(b[ind], Signal_anat[ind,0], label=r"$radius$ = " + str(round(radius[0],1)) +  " $\mu m$", color=colors[0])
ax2.plot(b[ind], Signal_anat[ind,1], label=r"$radius$ = " + str(round(radius[1],1)) +  " $\mu m$", color=colors[1])
ax2.plot(b[ind], Signal_anat[ind,2], label=r"$radius$ = " + str(round(radius[2],1)) +  " $\mu m$", color=colors[2])
ax2.set_yscale('log')
ax2.set_ylim([0.1, 1.1])
plt.xlabel(r"$b$" + " $(ms/ \mu m^2)$", fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

plt.savefig('Fig3_Diffraction_peaks_Delta_' + str(round(BigDelta6,2)) + '_delta_' + str(round(smalldelta6,2)) + '_D_'  + str(D) + '.png', bbox_inches='tight', dpi=600)
plt.show()
