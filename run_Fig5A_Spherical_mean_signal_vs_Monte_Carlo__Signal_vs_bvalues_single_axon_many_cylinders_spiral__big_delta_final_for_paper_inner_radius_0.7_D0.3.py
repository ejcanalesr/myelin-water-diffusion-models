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

b0      = bvalues == 0

# ---------------------- Define experimental parameters ------------------------
D          = 0.3 # um^2/ms
Nterms     = 40 # For the series
# Lori approximation
Lori        = True
pulse       = 'trapezoid'
#pulse       = 'rectangle'

# ---------------------------- Load MC simulations -----------------------------
g_ratio = 0.7
r_i     = 0.7
r_o     = r_i/g_ratio
radius_MC  = np.round(np.linspace(r_i, r_o, int(np.round((r_o - r_i)/0.1)+1) ), 1)

radius_discrete = np.round(radius_MC,1)
MC = np.zeros((6, radius_discrete.shape[0]))
for i in range(radius_discrete.shape[0]):
    radius_i  = radius_discrete[i]

    Signal_MC_protocol  = np.loadtxt('Data/MC_simulations_Connectome2_waveforms_6bvals/simulationsOuputsD03_G500_WF/myelin_' + str(radius_i) + '0_um_D_0.3_2D_G500_WF_DWI.txt')

    S0_protocol     = np.mean(Signal_MC_protocol[b0])
    SMT_protocol_b1 = np.mean(Signal_MC_protocol[ind_b1])/S0_protocol
    SMT_protocol_b2 = np.mean(Signal_MC_protocol[ind_b2])/S0_protocol
    SMT_protocol_b3 = np.mean(Signal_MC_protocol[ind_b3])/S0_protocol
    SMT_protocol_b4 = np.mean(Signal_MC_protocol[ind_b4])/S0_protocol
    SMT_protocol_b5 = np.mean(Signal_MC_protocol[ind_b5])/S0_protocol
    SMT_protocol_b6 = np.mean(Signal_MC_protocol[ind_b6])/S0_protocol

    Signal_MC = np.array([SMT_protocol_b1, SMT_protocol_b2, SMT_protocol_b3, SMT_protocol_b4, SMT_protocol_b5, SMT_protocol_b6])
    MC[:,i] = Signal_MC
    # --------------------------------------------------------------------------
#end

S_MC_06 = MC[:,0]
S_MC_07 = MC[:,1]
S_MC_08 = MC[:,2]
S_MC_09 = MC[:,3]

# ---------------------------- Load MC simulations -----------------------------
# --------------------------------- SPIRAL ----------------------------------- #
Signal_spiral = np.loadtxt('Data/MC_simulations_Connectome2_waveforms_6bvals/Spiral/simulations_ouputs_myelin_spiral_model_D03_WF/myelin_0.70_um_D_0.3_2D_spiral_wraps_rep_00_DWI.txt')

S0_spiral     = np.mean(Signal_spiral[b0])
SMT_spiral_b1 = np.mean(Signal_spiral[ind_b1])/S0_spiral
SMT_spiral_b2 = np.mean(Signal_spiral[ind_b2])/S0_spiral
SMT_spiral_b3 = np.mean(Signal_spiral[ind_b3])/S0_spiral
SMT_spiral_b4 = np.mean(Signal_spiral[ind_b4])/S0_spiral
SMT_spiral_b5 = np.mean(Signal_spiral[ind_b5])/S0_spiral
SMT_spiral_b6 = np.mean(Signal_spiral[ind_b6])/S0_spiral

# -------------------------Generate Signals -----------------------------------#
b1          = 0.80034088
Signal_radial1 = np.zeros((bvecs.shape[0], radius_MC.shape[0]))
for i in range(0, bvecs.shape[0]):
    # cos_alpha = [0 0 1]*bvecs[i,:]
    cos_alpha2  = bvecs[i,2]**2
    sin_alpha   = np.sqrt(1 -cos_alpha2)
    Signal_radial1[i,:] = Radial_signal(Nterms, D, b1, radius_MC, BigDelta1, smalldelta1, E1, sin_alpha, Lori, pulse) * np.exp(-b1 * D * cos_alpha2)
#end
Signal1 = np.mean(Signal_radial1, axis=0)

Signal1_anat   = SMT_signal(Nterms, D, b1, radius_MC, BigDelta1, smalldelta1, E1, Lori, pulse)
Signal1_anat_G = SMT_signal_Gaussian(D, b1, radius_MC, BigDelta1, smalldelta1, E1, Lori, pulse)

b2 =  0.99967705
Signal_radial2 = np.zeros((bvecs.shape[0], radius_MC.shape[0]))
for i in range(0, bvecs.shape[0]):
    cos_alpha2  = bvecs[i,2]**2
    sin_alpha   = np.sqrt(1 -cos_alpha2)
    Signal_radial2[i,:] = Radial_signal(Nterms, D, b2, radius_MC, BigDelta2, smalldelta2, E2, sin_alpha, Lori, pulse) * np.exp(-b2 * D * cos_alpha2)
#end
Signal2 = np.mean(Signal_radial2, axis=0)

Signal2_anat   = SMT_signal(Nterms, D, b2, radius_MC, BigDelta2, smalldelta2, E2, Lori, pulse)
Signal2_anat_G = SMT_signal_Gaussian(D, b2, radius_MC, BigDelta2, smalldelta2, E2, Lori, pulse)

b3 = 1.4996794
Signal_radial3 = np.zeros((bvecs.shape[0], radius_MC.shape[0]))
for i in range(0, bvecs.shape[0]):
    cos_alpha2  = bvecs[i,2]**2
    sin_alpha   = np.sqrt(1 -cos_alpha2)
    Signal_radial3[i,:] = Radial_signal(Nterms, D, b3, radius_MC, BigDelta3, smalldelta3, E3, sin_alpha, Lori, pulse) * np.exp(-b3 * D * cos_alpha2)
#end
Signal3 = np.mean(Signal_radial3, axis=0)

Signal3_anat   = SMT_signal(Nterms, D, b3, radius_MC, BigDelta3, smalldelta3, E3, Lori, pulse)
Signal3_anat_G = SMT_signal_Gaussian(D, b3, radius_MC, BigDelta3, smalldelta3, E3, Lori, pulse)

b4 = 2.0003624
Signal_radial4 = np.zeros((bvecs.shape[0], radius_MC.shape[0]))
for i in range(0, bvecs.shape[0]):
    cos_alpha2  = bvecs[i,2]**2
    sin_alpha   = np.sqrt(1 -cos_alpha2)
    Signal_radial4[i,:] = Radial_signal(Nterms, D, b4, radius_MC, BigDelta4, smalldelta4, E4, sin_alpha, Lori, pulse) * np.exp(-b4 * D * cos_alpha2)
#end
Signal4 = np.mean(Signal_radial4, axis=0)

Signal4_anat   = SMT_signal(Nterms, D, b4, radius_MC, BigDelta4, smalldelta4, E4, Lori, pulse)
Signal4_anat_G = SMT_signal_Gaussian(D, b4, radius_MC, BigDelta4, smalldelta4, E4, Lori, pulse)

b5 = 2.4995226
Signal_radial5 = np.zeros((bvecs.shape[0], radius_MC.shape[0]))
for i in range(0, bvecs.shape[0]):
    cos_alpha2  = bvecs[i,2]**2
    sin_alpha   = np.sqrt(1 -cos_alpha2)
    Signal_radial5[i,:] = Radial_signal(Nterms, D, b5, radius_MC, BigDelta5, smalldelta5, E5, sin_alpha, Lori, pulse) * np.exp(-b5 * D * cos_alpha2)
#end
Signal5 = np.mean(Signal_radial5, axis=0)

Signal5_anat   = SMT_signal(Nterms, D, b5, radius_MC, BigDelta5, smalldelta5, E5, Lori, pulse)
Signal5_anat_G = SMT_signal_Gaussian(D, b5, radius_MC, BigDelta5, smalldelta5, E5, Lori, pulse)

b6 = 3.0007138
Signal_radial6 = np.zeros((bvecs.shape[0], radius_MC.shape[0]))
for i in range(0, bvecs.shape[0]):
    cos_alpha2  = bvecs[i,2]**2
    sin_alpha   = np.sqrt(1 -cos_alpha2)
    Signal_radial6[i,:] = Radial_signal(Nterms, D, b6, radius_MC, BigDelta6, smalldelta6, E6, sin_alpha, Lori, pulse) * np.exp(-b6 * D * cos_alpha2)
#end
Signal6 = np.mean(Signal_radial6, axis=0)

Signal6_anat   = SMT_signal(Nterms, D, b6, radius_MC, BigDelta6, smalldelta6, E6, Lori, pulse)
Signal6_anat_G = SMT_signal_Gaussian(D, b6, radius_MC, BigDelta6, smalldelta6, E6, Lori, pulse)

# ------------------------ Plot results ----------------------------------------
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

plt.figure(figsize=(7.5,5))

b_values = np.array([b1, b2, b3, b4, b5, b6])

Signal_06_anat = np.array([Signal1_anat[0], Signal2_anat[0], Signal3_anat[0], Signal4_anat[0], Signal5_anat[0], Signal6_anat[0]])
Signal_07_anat = np.array([Signal1_anat[1], Signal2_anat[1], Signal3_anat[1], Signal4_anat[1], Signal5_anat[1], Signal6_anat[1] ])
Signal_08_anat = np.array([Signal1_anat[2], Signal2_anat[2], Signal3_anat[2], Signal4_anat[2], Signal5_anat[2], Signal6_anat[2]])
Signal_09_anat = np.array([Signal1_anat[3], Signal2_anat[3], Signal3_anat[3], Signal4_anat[3], Signal5_anat[3], Signal6_anat[3]])

Signal_06_anat_G = np.array([Signal1_anat_G[0], Signal2_anat_G[0], Signal3_anat_G[0], Signal4_anat_G[0], Signal5_anat_G[0], Signal6_anat_G[0]])
Signal_07_anat_G = np.array([Signal1_anat_G[1], Signal2_anat_G[1], Signal3_anat_G[1], Signal4_anat_G[1], Signal5_anat_G[1], Signal6_anat_G[1]])
Signal_08_anat_G = np.array([Signal1_anat_G[2], Signal2_anat_G[2], Signal3_anat_G[2], Signal4_anat_G[2], Signal5_anat_G[2], Signal6_anat_G[2]])
Signal_09_anat_G = np.array([Signal1_anat_G[3], Signal2_anat_G[3], Signal3_anat_G[3], Signal4_anat_G[3], Signal5_anat_G[3], Signal6_anat_G[3]])

S_MC_b1 = np.sum(MC[0,:] * radius_MC)/np.sum(radius_MC)
S_MC_b2 = np.sum(MC[1,:] * radius_MC)/np.sum(radius_MC)
S_MC_b3 = np.sum(MC[2,:] * radius_MC)/np.sum(radius_MC)
S_MC_b4 = np.sum(MC[3,:] * radius_MC)/np.sum(radius_MC)
S_MC_b5 = np.sum(MC[4,:] * radius_MC)/np.sum(radius_MC)
S_MC_b6 = np.sum(MC[5,:] * radius_MC)/np.sum(radius_MC)

S_MC_averague = np.array ([S_MC_b1, S_MC_b2, S_MC_b3, S_MC_b4, S_MC_b5, S_MC_b6])
S_MC_spiral   = np.array ([SMT_spiral_b1, SMT_spiral_b2, SMT_spiral_b3, SMT_spiral_b4, SMT_spiral_b5, SMT_spiral_b6])

m1, = plt.plot(b_values, Signal_06_anat, label=r"$radius$ = " + str(round(radius_MC[0],1)) +  " $\mu m$: Analytical (Eq.[23])", color=colors[0], linewidth=1)
m2, = plt.plot(b_values, Signal_07_anat, label=r"$radius$ = " + str(round(radius_MC[1],1)) +  " $\mu m$: Analytical (Eq.[23])", color=colors[1], linewidth=1)
m3, = plt.plot(b_values, Signal_08_anat, label=r"$radius$ = " + str(round(radius_MC[2],1)) +  " $\mu m$: Analytical (Eq.[23])", color=colors[2], linewidth=1)
m4, = plt.plot(b_values, Signal_09_anat, label=r"$radius$ = " + str(round(radius_MC[3],1)) +  " $\mu m$: Analytical (Eq.[23])", color=colors[3], linewidth=1)

g1, = plt.plot(b_values, Signal_06_anat_G, label=r"Gaussian (Eq.[24])", color=colors[0], linestyle='-.', linewidth=1)
g2, = plt.plot(b_values, Signal_07_anat_G, label=r"Gaussian (Eq.[24])", color=colors[1], linestyle='-.', linewidth=1)
g3, = plt.plot(b_values, Signal_08_anat_G, label=r"Gaussian (Eq.[24])", color=colors[2], linestyle='-.', linewidth=1)
g4, = plt.plot(b_values, Signal_09_anat_G, label=r"Gaussian (Eq.[24])", color=colors[3], linestyle='-.', linewidth=1)

mc1, = plt.plot(b_values, S_MC_06, label="MC signal", marker='o', color=colors[0], linewidth=0)
mc2, = plt.plot(b_values, S_MC_07, label="MC signal", marker='o', color=colors[1], linewidth=0)
mc3, = plt.plot(b_values, S_MC_08, label="MC signal", marker='o', color=colors[2], linewidth=0)
mc4, = plt.plot(b_values, S_MC_09, label="MC signal", marker='o', color=colors[3], linewidth=0)

spiral,     = plt.plot(b_values, S_MC_spiral, marker='x', color='black', mfc='none', linewidth=0, label=r"MC signal: Spiral geometry")
concentric, = plt.plot(b_values, S_MC_averague, marker='s', color='cyan', mfc='none', linewidth=0, label=r"MC signal: Multiple concentric cylinders (Eq.[19])")

ax = plt.gca()
ax.set_xlim([np.min(b_values) - 0.1, np.max(b_values) + 0.1])

ax.set_ylim([0.2, 1.0])
#ax.set_yscale('log')

leg1 = ax.legend(handles=[m1, m2 ,m3, m4, g1, g2, g3, g4, mc1, mc2, mc3, mc4],fontsize=10, ncol=3, handleheight=0, labelspacing=0.75, frameon=False, title=r"$D$ = " + str(D) + "$\mu m^2/ms$", title_fontsize=13, bbox_to_anchor=(0.5, 0.125), loc='lower center')
ax.add_artist(leg1)

ax.legend(handles=[spiral, concentric],fontsize=10, ncol=1, handleheight=0, labelspacing=0.75, frameon=False, bbox_to_anchor=(0.5, 0.0), loc='lower center')

plt.xticks(b_values, fontsize=11.5)
plt.yticks(fontsize=11.5)
#plt.yscale('log')
plt.ylabel(r"Normalized Signal Amplitude", fontsize=13)
plt.xlabel(r"$b-value$" + " $(ms/ \mu m^2)$", fontsize=14)

plt.title(r"Myelin Water Spherical Mean dMRI Signal", fontsize=14)
plt.savefig('Fig5A_Radius_SMT_Anat_MC_Delta_single_axon_multiple_layers' + '_D_'  + str(D) + '_inner_radius' + str(r_i) + '.png', bbox_inches='tight', dpi=300)

plt.show()
