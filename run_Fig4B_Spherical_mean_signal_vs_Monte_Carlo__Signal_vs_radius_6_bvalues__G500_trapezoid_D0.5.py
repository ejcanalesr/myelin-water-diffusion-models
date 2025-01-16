#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Erick Jorge Canales-Rodriguez, 2024

from myelin_water_diffusionMRI_models import *

# Load bvecs
hardi_fname, hardi_bval_fname, hardi_bvec_fname = get_fnames('stanford_hardi')
bvals, bvecs = read_bvals_bvecs(hardi_bval_fname, hardi_bvec_fname)
bvecs = bvecs[bvals>0,:]

# ---------------------- Define experimental parameters ------------------------
D          = 0.5 # um^2/ms
Nterms     = 40 # For the series
# Lori approximation
Lori        = True
pulse       = 'trapezoid'
#pulse       = 'rectangle'

# ---------------------------- Load MC simulations -----------------------------
radius  = np.linspace(0.1, 5.0, 50)
radius_discrete = np.round(radius,1)


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

MC = np.zeros((6, radius_discrete.shape[0]))
for i in range(radius_discrete.shape[0]):
    radius_i  = radius_discrete[i]
    Signal_MC_protocol  = np.loadtxt('Data/MC_simulations_Connectome2_waveforms_6bvals/simulationsOuputsD05_G500_WF/myelin_' + str(radius_i) + '0_um_D_0.5_2D_G500_WF_DWI.txt')

    S0_protocol     = np.mean(Signal_MC_protocol[b0])
    SMT_protocol_b1 = np.mean(Signal_MC_protocol[ind_b1])/S0_protocol
    SMT_protocol_b2 = np.mean(Signal_MC_protocol[ind_b2])/S0_protocol
    SMT_protocol_b3 = np.mean(Signal_MC_protocol[ind_b3])/S0_protocol
    SMT_protocol_b4 = np.mean(Signal_MC_protocol[ind_b4])/S0_protocol
    SMT_protocol_b5 = np.mean(Signal_MC_protocol[ind_b5])/S0_protocol
    SMT_protocol_b6 = np.mean(Signal_MC_protocol[ind_b6])/S0_protocol

    Signal_MC = np.array([SMT_protocol_b1, SMT_protocol_b2, SMT_protocol_b4, SMT_protocol_b3, SMT_protocol_b5, SMT_protocol_b6])
    MC[:,i] = Signal_MC
    # --------------------------------------------------------------------------
#end

S_MC_b1 = MC[0,:]
S_MC_b2 = MC[1,:]
S_MC_b4 = MC[2,:]
S_MC_b3 = MC[3,:]
S_MC_b5 = MC[4,:]
S_MC_b6 = MC[5,:]

# -------------------------Generate Signals -----------------------------------#
radius_cont     = np.linspace(0.001, 5, 50)

#b1   = 0.8
b1    = 0.80034088
Signal_radial1 = np.zeros((bvecs.shape[0], radius_cont.shape[0]))
for i in range(0, bvecs.shape[0]):
    # cos_alpha = [0 0 1]*bvecs[i,:]
    cos_alpha2  = bvecs[i,2]**2
    sin_alpha   = np.sqrt(1 -cos_alpha2)
    Signal_radial1[i,:] = Radial_signal(Nterms, D, b1, radius_cont, BigDelta1, smalldelta1, E1, sin_alpha, Lori, pulse) * np.exp(-b1 * D * cos_alpha2)
#end
Signal1 = np.mean(Signal_radial1, axis=0)

Signal1_anat   = SMT_signal(Nterms, D, b1, radius_cont, BigDelta1, smalldelta1, E1, Lori, pulse)
Signal1_anat_G = SMT_signal_Gaussian(D, b1, radius_cont, BigDelta1, smalldelta1, E1, Lori, pulse)

#b2 = 1.0
b2 =  0.99967705
Signal_radial2 = np.zeros((bvecs.shape[0], radius_cont.shape[0]))
for i in range(0, bvecs.shape[0]):
    cos_alpha2  = bvecs[i,2]**2
    sin_alpha   = np.sqrt(1 -cos_alpha2)
    Signal_radial2[i,:] = Radial_signal(Nterms, D, b2, radius_cont, BigDelta2, smalldelta2, E2, sin_alpha, Lori, pulse) * np.exp(-b2 * D * cos_alpha2)
#end
Signal2 = np.mean(Signal_radial2, axis=0)

Signal2_anat   = SMT_signal(Nterms, D, b2, radius_cont, BigDelta2, smalldelta2, E2, Lori, pulse)
Signal2_anat_G = SMT_signal_Gaussian(D, b2, radius_cont, BigDelta2, smalldelta2, E2, Lori, pulse)

#b3 = 1.5
b3 = 1.4996794
Signal_radial3 = np.zeros((bvecs.shape[0], radius_cont.shape[0]))
for i in range(0, bvecs.shape[0]):
    cos_alpha2  = bvecs[i,2]**2
    sin_alpha   = np.sqrt(1 -cos_alpha2)
    Signal_radial3[i,:] = Radial_signal(Nterms, D, b3, radius_cont, BigDelta3, smalldelta3, E3, sin_alpha, Lori, pulse) * np.exp(-b3 * D * cos_alpha2)
#end
Signal3 = np.mean(Signal_radial3, axis=0)

Signal3_anat   = SMT_signal(Nterms, D, b3, radius_cont, BigDelta3, smalldelta3, E3, Lori, pulse)
Signal3_anat_G = SMT_signal_Gaussian(D, b3, radius_cont, BigDelta3, smalldelta3, E3, Lori, pulse)

# b4  = 2.0
b4 = 2.0003624
Signal_radial4 = np.zeros((bvecs.shape[0], radius_cont.shape[0]))
for i in range(0, bvecs.shape[0]):
    cos_alpha2  = bvecs[i,2]**2
    sin_alpha   = np.sqrt(1 -cos_alpha2)
    Signal_radial4[i,:] = Radial_signal(Nterms, D, b4, radius_cont, BigDelta4, smalldelta4, E4, sin_alpha, Lori, pulse) * np.exp(-b4 * D * cos_alpha2)
#end
Signal4 = np.mean(Signal_radial4, axis=0)

Signal4_anat   = SMT_signal(Nterms, D, b4, radius_cont, BigDelta4, smalldelta4, E4, Lori, pulse)
Signal4_anat_G = SMT_signal_Gaussian(D, b4, radius_cont, BigDelta4, smalldelta4, E4, Lori, pulse)

#b5 = 2.5
b5 = 2.4995226
Signal_radial5 = np.zeros((bvecs.shape[0], radius_cont.shape[0]))
for i in range(0, bvecs.shape[0]):
    cos_alpha2  = bvecs[i,2]**2
    sin_alpha   = np.sqrt(1 -cos_alpha2)
    Signal_radial5[i,:] = Radial_signal(Nterms, D, b5, radius_cont, BigDelta5, smalldelta5, E5, sin_alpha, Lori, pulse) * np.exp(-b5 * D * cos_alpha2)
#end
Signal5 = np.mean(Signal_radial5, axis=0)

Signal5_anat   = SMT_signal(Nterms, D, b5, radius_cont, BigDelta5, smalldelta5, E5, Lori, pulse)
Signal5_anat_G = SMT_signal_Gaussian(D, b5, radius_cont, BigDelta5, smalldelta5, E5, Lori, pulse)

#b6 = 3.0
b6 = 3.0007138
Signal_radial6 = np.zeros((bvecs.shape[0], radius_cont.shape[0]))
for i in range(0, bvecs.shape[0]):
    cos_alpha2  = bvecs[i,2]**2
    sin_alpha   = np.sqrt(1 -cos_alpha2)
    Signal_radial6[i,:] = Radial_signal(Nterms, D, b6, radius_cont, BigDelta6, smalldelta6, E6, sin_alpha, Lori, pulse) * np.exp(-b6 * D * cos_alpha2)
#end
Signal6 = np.mean(Signal_radial6, axis=0)

Signal6_anat   = SMT_signal(Nterms, D, b6, radius_cont, BigDelta6, smalldelta6, E6, Lori, pulse)
Signal6_anat_G = SMT_signal_Gaussian(D, b6, radius_cont, BigDelta6, smalldelta6, E6, Lori, pulse)

# ------------------------ Plot results ----------------------------------------
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

#plt.figure(figsize=(7.5,5))
plt.figure(figsize=(7.5,3.9))

legend = False
if legend == False:
    plt.plot(radius_cont, Signal1_anat, color=colors[0])
    plt.plot(radius_cont, Signal2_anat, color=colors[1])
    plt.plot(radius_cont, Signal3_anat, color=colors[2])
    plt.plot(radius_cont, Signal4_anat, color=colors[3])
    plt.plot(radius_cont, Signal5_anat, color=colors[4])
    plt.plot(radius_cont, Signal6_anat, color=colors[5])

    plt.plot(radius_cont, Signal1_anat_G, color=colors[0], linestyle='-.')
    plt.plot(radius_cont, Signal2_anat_G, color=colors[1], linestyle='-.')
    plt.plot(radius_cont, Signal3_anat_G, color=colors[2], linestyle='-.')
    plt.plot(radius_cont, Signal4_anat_G, color=colors[3], linestyle='-.')
    plt.plot(radius_cont, Signal5_anat_G, color=colors[4], linestyle='-.')
    plt.plot(radius_cont, Signal6_anat_G, color=colors[5], linestyle='-.')

    plt.plot(radius_discrete, S_MC_b1, marker='o', color=colors[0], linewidth=0, markersize=4)
    plt.plot(radius_discrete, S_MC_b2, marker='o', color=colors[1], linewidth=0, markersize=4)
    plt.plot(radius_discrete, S_MC_b3, marker='o', color=colors[2], linewidth=0, markersize=4)
    plt.plot(radius_discrete, S_MC_b4, marker='o', color=colors[3], linewidth=0, markersize=4)
    plt.plot(radius_discrete, S_MC_b5, marker='o', color=colors[4], linewidth=0, markersize=4)
    plt.plot(radius_discrete, S_MC_b6, marker='o', color=colors[5], linewidth=0, markersize=4)
elif legend == True:
    plt.plot(radius_cont, Signal1_anat, label=r"$b$ = " + str(round(b1,1)) +  " $ms/ \mu m^2$: Analytical (Eq.[23])", color=colors[0])
    plt.plot(radius_cont, Signal2_anat, label=r"$b$ = " + str(round(b2,1)) +  " $ms/ \mu m^2$: Analytical (Eq.[23])", color=colors[1])
    plt.plot(radius_cont, Signal3_anat, label=r"$b$ = " + str(round(b3,1)) +  " $ms/ \mu m^2$: Analytical (Eq.[23])", color=colors[2])
    plt.plot(radius_cont, Signal4_anat, label=r"$b$ = " + str(round(b4,1)) +  " $ms/ \mu m^2$: Analytical (Eq.[23])", color=colors[3])
    plt.plot(radius_cont, Signal5_anat, label=r"$b$ = " + str(round(b5,1)) +  " $ms/ \mu m^2$: Analytical (Eq.[23])", color=colors[4])
    plt.plot(radius_cont, Signal6_anat, label=r"$b$ = " + str(round(b6,1)) +  " $ms/ \mu m^2$: Analytical (Eq.[23])", color=colors[5])

    plt.plot(radius_cont, Signal1_anat_G, label=r"Gaussian (Eq.[24])", color=colors[0], linestyle='-.')
    plt.plot(radius_cont, Signal2_anat_G, label=r"Gaussian (Eq.[24])", color=colors[1], linestyle='-.')
    plt.plot(radius_cont, Signal3_anat_G, label=r"Gaussian (Eq.[24])", color=colors[2], linestyle='-.')
    plt.plot(radius_cont, Signal4_anat_G, label=r"Gaussian (Eq.[24])", color=colors[3], linestyle='-.')
    plt.plot(radius_cont, Signal5_anat_G, label=r"Gaussian (Eq.[24])", color=colors[4], linestyle='-.')
    plt.plot(radius_cont, Signal6_anat_G, label=r"Gaussian (Eq.[24])", color=colors[5], linestyle='-.')

    plt.plot(radius_discrete, S_MC_b1, label="MC signal", marker='o', color=colors[0], linewidth=0, markersize=4)
    plt.plot(radius_discrete, S_MC_b2, label="MC signal", marker='o', color=colors[1], linewidth=0, markersize=4)
    plt.plot(radius_discrete, S_MC_b3, label="MC signal", marker='o', color=colors[2], linewidth=0, markersize=4)
    plt.plot(radius_discrete, S_MC_b4, label="MC signal", marker='o', color=colors[3], linewidth=0, markersize=4)
    plt.plot(radius_discrete, S_MC_b5, label="MC signal", marker='o', color=colors[4], linewidth=0, markersize=4)
    plt.plot(radius_discrete, S_MC_b6, label="MC signal", marker='o', color=colors[5], linewidth=0, markersize=4)
#end

ax = plt.gca()
#ax.set_xlim([0, 1.5])
ax.set_ylim([0.0, 1])
ax.legend(fontsize=10, ncol=3, handleheight=1.5, labelspacing=0.5, frameon=False, title=r"$D$ = " + str(D) + " $\mu m^2/ms$", title_fontsize=12)
plt.xticks(fontsize=11.5)
plt.yticks(fontsize=11.5)
#plt.yscale('log')
plt.ylabel(r"Normalized Signal Amplitude", fontsize=13)
plt.xlabel(r"$radius$" + " $(\mu m)$", fontsize=14)

plt.xticks(np.linspace(0,5,11))

plt.title(r"Myelin Water Spherical Mean dMRI Signal", fontsize=14)
plt.savefig('Fig4B_Radius_Delta_SMT_Anat_MC_Delta' + '_D_'  + str(D) +  '_G500_' + pulse + '.png', bbox_inches='tight', dpi=300)
plt.show()
