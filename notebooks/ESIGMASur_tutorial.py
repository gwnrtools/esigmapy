#!/usr/bin/env python
# coding: utf-8

# In[8]:


get_ipython().system('pip install ~/src/esigmapy/')


# # ESIGMASur tutorial

# This notebook illustrates the usage of `ESIGMASur`, which is a *non-spinning* *time-domain* surrogate model of the *(2,2)-mode* of `ESIGMAHM`.

# > **Note:** If you are evaluating this notebook in its original location within the repository, the cell below will correctly set the path of the surrogate data files to the `bash` environment variable `ESIGMASUR_DATA_PATH`. However, if you have moved this notebook, please point `ESIGMASUR_DATA_PATH` to the path of the [surrogate data directory](https://github.com/gwnrtools/esigmapy/tree/master/esigmapy/surrogate/data) on your local machine. 

# In[1]:


import os

# Required for efficient surrogate evaluation by avoiding multi-threading overheads
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import matplotlib.pyplot as plt
import esigmapy.surrogate as esigmasur

# Configuring some plot settings
plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "font.size": 12,
    }
)

# Point bash variable ESIGMASUR_DATA_PATH to the directory where surrogate data is stored (esigmapy/surrogate/data)
sur_data_dir = "../esigmapy/surrogate/data"
sur_data_dir = os.path.abspath(sur_data_dir)
os.environ["ESIGMASUR_DATA_PATH"] = sur_data_dir


# ## 1. InspiralESIGMASur
# 
# `InspiralESIGMASur` is a *non-spinning* *time-domain* surrogate model of `InspiralESIGMA`, the *inspiral* *(2,2)-mode* piece of `ESIGMAHM`. The surrogate is built using a new scalable eccentric surrogate modeling technique presented in [arxiv:2510.00116](https://arxiv.org/abs/2510.00116).
# 
# `InspiralESIGMASur` can generate eccentric inspiral waveforms of lengths up to $2.77 \times 10^{6}M$ for binaries with mass-ratios $q:=m_1/m_2 \in [1,6]$, and reference eccentricities and mean anomalies $e_{\rm{ref}} \in [0, 0.431]$ and $l_{\rm{ref}} \in [0, 2 \pi)$ measured at a reference time $t_{\rm{ref}} = -2.75 \times 10^6M$ before the end of the inspiral ($t=0$).

# ### 1.1 Waveform polarizations
# The polarizations $h_+$ and $h_\times$ can be generated via the `get_inspiral_esigmasur_waveform` function. They are returned as `PyCBC` `TimeSeries` objects.

# In[2]:


m1 = 7.0  # masses (in solar masses)
m2 = 3.0
reference_eccentricity = 0.43  # reference eccentricity
reference_mean_anomaly = 60 * np.pi / 180.0  # reference mean anomaly

distance = 400.0  # source luminosity distance (in Mpc)
inclination = 30 * np.pi / 180.0  # orbital inclination with line-of-sight

delta_t = 1 / 2**12  # time grid-spacing (in s)
# Waveform starting time (in s); 2.77e6M (the maximum surrogate length) corresponds
# to roughly 136s for a 10 Msun binary. t=0 corresponds to the end of waveform,
# so the starting time should be a negative real number
t_start = -136.0

hp, hc = esigmasur.get_inspiral_esigmasur_waveform(
    mass1=m1,
    mass2=m2,
    reference_eccentricity=reference_eccentricity,
    reference_mean_anomaly=reference_mean_anomaly,
    delta_t=delta_t,
    t_start=t_start,
    distance=distance,
    inclination=inclination,
)

plt.figure(figsize=(10, 4))
plt.title(
    rf"""$m_1={m1} M_\odot, m_2={m2} M_\odot, e_{{\rm{{ref}}}}={reference_eccentricity}, l_{{\rm{{ref}}}}={reference_mean_anomaly:.2f}, t_{{\rm{{start}}}}={t_start:.2f}\rm{{s}}$
          $d_L={distance}\rm{{Mpc}}, \iota={inclination:.2f}$"""
)
hp.plot(label=r"$h_+$")
hc.plot(label=r"$h_\times$")
plt.xlabel(r"$t (s)$")
plt.ylabel(r"$h$")
plt.legend()
plt.tight_layout()


# In[2]:


m1 = 7.0  # masses (in solar masses)
m2 = 3.0
reference_eccentricity = 0.43  # reference eccentricity
reference_mean_anomaly = 60 * np.pi / 180.0  # reference mean anomaly

distance = 400.0  # source luminosity distance (in Mpc)
inclination = 30 * np.pi / 180.0  # orbital inclination with line-of-sight

delta_t = 1 / 2**12  # time grid-spacing (in s)
# Waveform starting time (in s); 2.77e6M (the maximum surrogate length) corresponds
# to roughly 136s for a 10 Msun binary. t=0 corresponds to the end of waveform,
# so the starting time should be a negative real number
t_start = -136.0

hp, hc = esigmasur.get_inspiral_esigmasur_waveform(
    mass1=m1,
    mass2=m2,
    reference_eccentricity=reference_eccentricity,
    reference_mean_anomaly=reference_mean_anomaly,
    delta_t=delta_t,
    t_start=t_start,
    distance=distance,
    inclination=inclination,
)

plt.figure(figsize=(10, 4))
plt.title(
    rf"""$m_1={m1} M_\odot, m_2={m2} M_\odot, e_{{\rm{{ref}}}}={reference_eccentricity}, l_{{\rm{{ref}}}}={reference_mean_anomaly:.2f}, t_{{\rm{{start}}}}={t_start:.2f}\rm{{s}}$
          $d_L={distance}\rm{{Mpc}}, \iota={inclination:.2f}$"""
)
hp.plot(label=r"$h_+$")
hc.plot(label=r"$h_\times$")
plt.xlabel(r"$t (s)$")
plt.ylabel(r"$h$")
plt.legend()
plt.tight_layout()


# The first evaluation will be slower due to two one-time costs: (1) loading surrogate data from disk, and (2) `numba` JIT compilation of some routines.

# ### 1.2 Waveform modes
# One can also generate the spin-weighted spherical harmonic modes via the `get_inspiral_esigmasur_modes` function. Only the $(2, \pm2)$ modes are supported currently.

# In[3]:


modes = esigmasur.get_inspiral_esigmasur_modes(
    mass1=m1,
    mass2=m2,
    reference_eccentricity=reference_eccentricity,
    reference_mean_anomaly=reference_mean_anomaly,
    delta_t=delta_t,
    t_start=t_start,
    distance=distance,
)

plt.figure(figsize=(10, 4))
plt.title(
    rf"""$m_1={m1} M_\odot, m_2={m2} M_\odot, e_{{\rm{{ref}}}}={reference_eccentricity}, l_{{\rm{{ref}}}}={reference_mean_anomaly:.2f}, t_{{\rm{{start}}}}={t_start:.2f}\rm{{s}}$
          $d_L={distance}\rm{{Mpc}}$"""
)

mode_name = (2, 2)
ell, m = mode_name
plt.plot(
    modes[mode_name].sample_times.data,
    modes[mode_name].real().data,
    label=rf"$\Re(h_{{{ell} {m}}})$",
)
plt.plot(
    modes[mode_name].sample_times.data,
    modes[mode_name].imag().data,
    label=rf"$\Im(h_{{{ell} {m}}})$",
)
plt.legend(loc=2)
plt.xlabel(r"$t (s)$")
plt.ylabel(r"$h$")
plt.tight_layout()


# In[7]:


modes = esigmasur.get_inspiral_esigmasur_modes(
    mass1=m1,
    mass2=m2,
    reference_eccentricity=reference_eccentricity,
    reference_mean_anomaly=reference_mean_anomaly,
    delta_t=delta_t,
    t_start=t_start,
    distance=distance,
)

plt.figure(figsize=(10, 4))
plt.title(
    rf"""$m_1={m1} M_\odot, m_2={m2} M_\odot, e_{{\rm{{ref}}}}={reference_eccentricity}, l_{{\rm{{ref}}}}={reference_mean_anomaly:.2f}, t_{{\rm{{start}}}}={t_start:.2f}\rm{{s}}$
          $d_L={distance}\rm{{Mpc}}$"""
)

mode_name = (2, 2)
ell, m = mode_name
plt.plot(
    modes[mode_name].sample_times.data,
    modes[mode_name].real().data,
    label=rf"$\Re(h_{{{ell} {m}}})$",
)
plt.plot(
    modes[mode_name].sample_times.data,
    modes[mode_name].imag().data,
    label=rf"$\Im(h_{{{ell} {m}}})$",
)
plt.legend(loc=2)
plt.xlabel(r"$t (s)$")
plt.ylabel(r"$h$")
plt.tight_layout()


# ### 1.3 Evolution of orbital elements
# The surrogates of binary's orbital elements' evolution can also be accessed via the argument `return_orbital_params` in all of the above discussed waveform/mode functions. The available orbital elements are
# 
# - $e$: Orbital eccentricity
# - $l$: Mean anomaly
# - $x$: The post-Newtonian (PN) parameter. It's related to the orbit-averaged (azimuthal) orbital frequency
# 
# $e$ and $l$ surrogates are internally required for surrogate waveform generation to get the values of eccentricity and mean anomaly at empirical interpolation (EI) nodes (i.e., $e_{\rm{EI}}$ and $l_{\rm{EI}}$); see [arxiv:2510.00116](https://arxiv.org/abs/2510.00116) for details. $x$ surrogate is required by the merger-ringdown attachment algorithm (as detailed in [arxiv:2409.13866](https://arxiv.org/abs/2409.13866)) for producing an IMR waveform using the inspiral surrogate (demonstrated later in this notebook).

# In[8]:


orb_params_list = ["e", "l", "x"]
orb_vars, hp, hc = esigmasur.get_inspiral_esigmasur_waveform(
    mass1=m1,
    mass2=m2,
    reference_eccentricity=reference_eccentricity,
    reference_mean_anomaly=reference_mean_anomaly,
    delta_t=delta_t,
    t_start=t_start,
    distance=distance,
    inclination=inclination,
    return_orbital_params=orb_params_list,
)

fig, axs = plt.subplots(len(orb_params_list), sharex=True, figsize=(10, 10))
axs[0].set_title(
    rf"""$m_1={m1} M_\odot, m_2={m2} M_\odot, e_{{\rm{{ref}}}}={reference_eccentricity}, l_{{\rm{{ref}}}}={reference_mean_anomaly:.2f}, t_{{\rm{{start}}}}={t_start}\rm{{s}}$"""
)
for i, orb_params_name in enumerate(orb_params_list):
    axs[i].plot(
        orb_vars[orb_params_name].sample_times.data,
        orb_vars[orb_params_name].data,
        label=rf"{orb_params_name}",
    )
    axs[i].legend(loc=2)
plt.xlabel(r"$t (s)$")


# ---

# ## 2. IMRESIGMASur
# The inspiral surrogate `InspiralESIGMASur` can be used as a drop-in replacement for `InspiralESIGMA`, the eccentric inspiral piece in the [`ESIGMA` framework](https://arxiv.org/abs/2409.13866), and can be smoothly attached to a quasi-circular plunge-merger-ringdown piece (`NRSur7dq4`, by default) to produce a complete inspiral-merger-ringdown (IMR) model. We call this hybrid IMR surrogate model `IMRESIGMASur`. 
# 
# **Note:** This is not a single eccentric IMR surrogate, but rather a hybrid of two surrogates: an eccentric inspiral surrogate (`InspiralESIGMASur`) and a quasi-circular plunge-merger-ringdown surrogate (`NRSur7dq4`).   
# 
# > **Note:** Running the following two cells will require installing the surrogate data files for `NRSur7dq4` (see the instructions at [`ESIGMAPy`'s wiki](https://github.com/gwnrtools/esigmapy)).

# ### 2.1 Waveform polarizations

# In[9]:


hp, hc = esigmasur.get_imr_esigmasur_waveform(
    mass1=m1,
    mass2=m2,
    reference_eccentricity=reference_eccentricity,
    reference_mean_anomaly=reference_mean_anomaly,
    delta_t=delta_t,
    t_start=t_start,
    distance=distance,
    inclination=inclination,
)

plt.figure(figsize=(10, 4))
plt.title(
    rf"""$m_1={m1} M_\odot, m_2={m2} M_\odot, e_{{\rm{{ref}}}}={reference_eccentricity}, l_{{\rm{{ref}}}}={reference_mean_anomaly:.2f}, t_{{\rm{{start}}}}={t_start:.2f}\rm{{s}}$
          $d_L={distance}\rm{{Mpc}}, \iota={inclination:.2f}$"""
)
hp.plot(label=r"$h_+$")
hc.plot(label=r"$h_\times$")
plt.xlabel(r"$t (s)$")
plt.ylabel(r"$h$")
plt.legend()
plt.tight_layout()


# ### 2.2 Waveform Modes

# In[10]:


modes = esigmasur.get_imr_esigmasur_mode(
    mass1=m1,
    mass2=m2,
    reference_eccentricity=reference_eccentricity,
    reference_mean_anomaly=reference_mean_anomaly,
    delta_t=delta_t,
    t_start=t_start,
    distance=distance,
)

mode_name = (2, 2)
ell, m = mode_name
plt.figure(figsize=(10, 4))
plt.title(
    rf"""$m_1={m1} M_\odot, m_2={m2} M_\odot, e_{{\rm{{ref}}}}={reference_eccentricity}, l_{{\rm{{ref}}}}={reference_mean_anomaly:.2f}, t_{{\rm{{start}}}}={t_start}\rm{{s}}$
          $d_L={distance}\rm{{Mpc}}$"""
)
plt.plot(
    modes[mode_name].sample_times.data,
    modes[mode_name].real().data,
    label=rf"$\Re(h_{{{ell} {m}}})$",
)
plt.plot(
    modes[mode_name].sample_times.data,
    modes[mode_name].imag().data,
    label=rf"$\Im(h_{{{ell} {m}}})$",
)
plt.legend(loc=2)
plt.xlabel(r"$t (s)$")
plt.tight_layout()


# <hr style="border: 2px solid #555; margin: 20px 0;">

# ## 3. Advanced features and demonstrations

# ### 3.1 Evaluation over generic time-grids
# The inspiral surrogate `InspiralESIGMASur` can also be evaluated at a user-specified time-grid (sorted in ascending order), which can be *non-uniform*. This time-grid can be supplied to the `times` argument in the surrogate's waveform/mode generation functions as a `NumPy` array. 
# 
# As an example, we use this feature to sample the waveform uniformly in *GW phase* instead of time. Consequently, in the plot below one can observe that
# 1. within any cycle, the time sampling is denser near the peaks and troughs, where the waveform changes rapidly
# 2. the time-sampling becomes denser near the end of inspiral, where the waveform evolves rapidly

# In[11]:


# Generating the same waveform as above
hp, hc = esigmasur.get_inspiral_esigmasur_waveform(
    mass1=m1,
    mass2=m2,
    reference_eccentricity=reference_eccentricity,
    reference_mean_anomaly=reference_mean_anomaly,
    delta_t=delta_t,
    t_start=t_start,
    distance=distance,
    inclination=inclination,
)

# Computing a non-uniform time-grid that is uniformly sampled in GW phase instead
phase_gw = np.unwrap(np.arctan2(hc.data, hp.data))
num_points_per_cycle = 18
uniform_phase_grid = np.linspace(
    phase_gw[2],
    phase_gw[-2],
    int(num_points_per_cycle * (phase_gw[-2] - phase_gw[2]) / (2 * np.pi)),
    endpoint=True,
)
time_grid_uniformly_sampled_in_phase = np.interp(
    uniform_phase_grid, phase_gw, hp.sample_times.data
)

# Generating the waveform on the non-uniform time grid by passing it to the `times` argument
t, hp, hc = esigmasur.get_inspiral_esigmasur_waveform(
    mass1=m1,
    mass2=m2,
    reference_eccentricity=reference_eccentricity,
    reference_mean_anomaly=reference_mean_anomaly,
    times=time_grid_uniformly_sampled_in_phase,
    distance=distance,
    inclination=inclination,
    return_pycbc_timeseries=False,  # When giving a custom time-grid, can't return waveform as a PyCBC TimeSeries object
)

# Plotting
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
fig.suptitle(
    r"Waveform on a non-uniform time-grid that is uniformly sampled in GW phase"
)

ax1.plot(t, hp, marker=".", lw=0.9, markersize=5, label=r"$h_+$")
ax1.plot(t, hc, marker=".", lw=0.9, markersize=5, label=r"$h_\times$")
ax1.set_ylabel("$h$")
ax1.set_xlim(-136.0, -135.0)
ax1.legend()

ax2.plot(t, hp, marker=".", lw=0.9, markersize=5, label=r"$h_+$")
ax2.plot(t, hc, marker=".", lw=0.9, markersize=5, label=r"$h_\times$")
ax2.set_ylabel("$h$")
ax2.set_xlim(-0.1, 0.0)
ax2.legend()

plt.xlabel(r"$t (s)$")
plt.tight_layout()


# ### 3.2 Surrogate metrics
# One can explore the internals of `InspiralESIGMASur` by accessing the base surrogate object via the `get_surrogate_object` function.
# 
# As an example, we access the parameter space points at which the surrogate was trained from the surrogate's metadata file and print their total number. We also access the number of basis functions required by the various individual surrogate data-pieces (see [arxiv:2510.00116](https://arxiv.org/abs/2510.00116) for their meanings) by accessing their respective "B-matrices" produced by the Empirical Interpolation Method (EIM) (see Eq. 19 of [arXiv:1308.3565](https://arxiv.org/abs/1308.3565)).

# In[6]:


# Expose the base surrogate object
sur = esigmasur.get_surrogate_object()
# Surrogate's training parameter space: [q, e_ref, l_ref]
param_space = sur.get_metadata("training_param_space")

print(f"The surrogate is trained over {len(param_space)} InspiralESIGMA waveforms.\n")

full_eccentric_data_piece_names = {
    "res_amp": "Residual amplitude",
    "res_phase": "Residual phase",
    "res_circ_phase": "Residual circular phase",
    "e": "Eccentricity dynamics",
    "shifted_mean_anomaly": "Shifted mean anomaly dynamics",
}

full_circular_data_piece_names = {
    "amp": "Circular amplitude",
    "phase": "Circular phase",
}

print(f"The surrogate uses the following number of basis functions:")
for key, value in full_eccentric_data_piece_names.items():
    # Accessing the EIM-B matrices to find the number of basis functions
    print(f"{value}: {len(sur.eim_B[key])}")

for key, value in full_circular_data_piece_names.items():
    print(f"{value}: {len(sur.circ_sur.eim_B[key])}")


# ### 3.3 Surrogate's domain of validity in starting GW frequency
# `InspiralESIGMASur` can generate waveforms of lengths up to $2.77 \times 10^6M$ in time. However, in order to better understand the validity of the surrogate in terms of starting GW frequencies, we plot the minimum frequency with which these waveforms can be started for a $10M_\odot$ binary across the surrogate's training parameter space. Equivalently, one can also think in terms of the minimum binary mass for which the waveforms can be started from a fixed frequency, say $15\rm{Hz}$.
# 
# We find these metrics below by utilizing the starting values of the PN parameter $x$ stored in the surrogate's metadata at the training parameter space points.
# 
# **Note:** Since the frequencies are extracted from the PN parameter $x$, these correspond to the orbit-averaged (azimuthal) GW frequency, and NOT the instantaneous (2,2)-mode frequency.

# In[7]:


from esigmapy.utils import f22_from_x

# Total binary mass (in M_Sun) for which to compute the frequencies
M_eval = 10
# Values of the PN parameter "x" at the starting of the surrogate's full length
x_min_array = sur.get_metadata("x_min_array")
f_min_at_M_eval = f22_from_x(np.max(x_min_array), M=M_eval)
print(
    f"The minimum GW frequency till which the surrogate can generate waveforms for a {M_eval:.1f}M_Sun binary across its entire parameter space = {f_min_at_M_eval:.1f}Hz"
)

f_compute = 15.0
M_min_at_fcompute = M_eval * (f_min_at_M_eval / f_compute)
print(
    f"Equivalently, the minimum total binary mass for which the surrogate can generate waveforms from {f_compute}Hz across its entire parameter space = {M_min_at_fcompute:.1f}M_Sun"
)

# Converting PN parameter to the orbit-averaged GW frequency
plot_data = f22_from_x(x_min_array, M=M_eval)

fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(111, projection="3d")

img = ax.scatter(
    param_space[:, 1],
    param_space[:, 2],
    param_space[:, 0],
    c=plot_data,
    cmap="viridis",
    s=25,
    alpha=0.7,
)

cbar = fig.colorbar(img, ax=ax, shrink=0.7, pad=0.1)
cbar.set_label(r"$f_{\rm{start}} \, (\rm{Hz})$")

ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False

ax.set_title(rf"$f_{{\rm{{start}}}}$ @ ${M_eval:.0f} M_\odot$")
ax.set_zlabel(r"$q$")
ax.set_xlabel(r"$e_{\rm{ref}}$")
ax.set_ylabel(r"$l_{\rm{ref}}$")

plt.tight_layout()
plt.show()


# ### 3.4 Surrogate's evaluation speed  

# We compute the evaluation speed of `InspiralESIGMASur` at random points in mass-ratio $q = m_1/m_2$, reference eccentricity $e_{\rm{ref}}$ and reference mean anomaly $l_{\rm{ref}}$. The surrogate is started such that its starting frequency is $15\rm{Hz}$, and the waveforms are sampled at $4096 \rm{Hz}$.
# 
# > **Note:** The following cell may take about a minute to run.

# In[8]:


import time
from esigmapy.utils import x_from_f22

num_evals_per_mass = 25
# Starting frequency (in Hz) from which to generate the surrogate waveforms for timing.
f_start = 15.0
delta_t = 1 / 2**12
total_mass_array = np.asarray([20, 30, 40, 50, 60, 80, 100])

NLOOP = 10


def speeds(M, q, e_ref, l_ref, f_start, delta_t):
    m1 = q * M / (1 + q)
    m2 = M / (1 + q)
    t, orb_vars, modes = esigmasur.get_inspiral_esigmasur_modes(
        mass1=m1,
        mass2=m2,
        reference_eccentricity=e_ref,
        reference_mean_anomaly=l_ref,
        delta_t=delta_t,
        t_start=None,  # Generate the surrogate with its full duration
        return_pycbc_timeseries=False,
        return_orbital_params=["x"],
    )
    # Finding the time at which the surrogate reaches
    # the GW frequency of f_start
    x_start = x_from_f22(f_start, M=M)
    idx = np.argmax(orb_vars["x"] >= x_start)
    t_start = t[idx]

    t1 = time.perf_counter()
    t, modes = esigmasur.get_inspiral_esigmasur_modes(
        mass1=m1,
        mass2=m2,
        reference_eccentricity=e_ref,
        reference_mean_anomaly=l_ref,
        delta_t=delta_t,
        # Generate the surrogate with a duration such that
        # it starts from a GW frequency of f_start
        t_start=t_start,
        return_pycbc_timeseries=False,
    )
    sur_time = time.perf_counter() - t1

    if sur_time < 2.5e-2:
        temp_time = 0
        for _ in range(NLOOP - 1):
            t1 = time.perf_counter()
            t, modes = esigmasur.get_inspiral_esigmasur_modes(
                mass1=m1,
                mass2=m2,
                reference_eccentricity=e_ref,
                reference_mean_anomaly=l_ref,
                delta_t=delta_t,
                # Generate the surrogate with a duration such that
                # it starts from a GW frequency of f_start
                t_start=t_start,
                distance=distance,
                return_pycbc_timeseries=False,
            )
            temp_time += time.perf_counter() - t1
        sur_time += temp_time
        sur_time /= NLOOP
    return sur_time


total_mass_array = total_mass_array[total_mass_array > M_min_at_fcompute]
total_mass_array = np.r_[M_min_at_fcompute, total_mass_array]
sur_eval_times_dict = {}

for total_mass in total_mass_array:
    rng = np.random.default_rng(seed=37)
    q_array = rng.uniform(1, 6, size=num_evals_per_mass)
    e_array = rng.uniform(0, 0.43, size=num_evals_per_mass)
    l_array = rng.uniform(0, 2 * np.pi, size=num_evals_per_mass)

    sur_eval_times = []
    for i in range(num_evals_per_mass):
        sur_eval_times.append(
            speeds(
                M=total_mass,
                q=q_array[i],
                e_ref=e_array[i],
                l_ref=l_array[i],
                f_start=f_start,
                delta_t=delta_t,
            )
        )
    sur_eval_times_dict[total_mass] = sur_eval_times
    print(f"Timing completed for M = {total_mass:.2f} M_Sun")


# In[11]:


import matplotlib.ticker as mticker

fig, ax = plt.subplots()
ax.set_title(
    rf"Surrogate evaluation time for $f_{{\rm{{start}}}} = {f_start} \rm{{Hz}}$, $f_{{\rm{{s}}}} = {1/delta_t:.0f} \rm{{Hz}}$"
)
ax.set_xlabel("Total mass ($M_\\odot$)")
ax.set_ylabel("Evaluation time (ms)")
ax.set_yscale("log")

alpha_scatter = 0.6
alpha_fill = 0.2
color = "dodgerblue"

for mass, eval_times in sur_eval_times_dict.items():
    eval_times = np.asarray(eval_times)

sur_eval_times_array = np.asarray(list(sur_eval_times_dict.values()))
sur_eval_times_median = np.median(sur_eval_times_array, axis=-1)
sur_speed_min = np.min(sur_eval_times_array, axis=-1)
sur_speed_max = np.max(sur_eval_times_array, axis=-1)

ax.plot(
    total_mass_array,
    sur_eval_times_median * 1.0e3,
    ls="--",
    color=color,
    marker=".",
    alpha=alpha_scatter,
    label="Median evaluation time",
)
ax.fill_between(
    total_mass_array,
    sur_speed_min * 1.0e3,
    sur_speed_max * 1.0e3,
    color=color,
    alpha=alpha_fill,
    label="Minimum/maximum evaluation time",
)


# To display only the first tick with a decimal entry
def mixed_formatter(x, pos):
    # pos = index of tick in visible ticks
    if pos == 0:
        return f"{x:.1f}"  # first tick: 1 decimal place
    else:
        return f"{x:.0f}"  # others: no decimals


ax.xaxis.set_major_formatter(mticker.FuncFormatter(mixed_formatter))
plt.xticks(total_mass_array)
plt.legend(frameon=False)
plt.tight_layout()

import subprocess

cpu_name = (
    subprocess.check_output("lscpu | grep 'Model name'", shell=True)
    .decode()
    .split(":")[1]
    .strip()
)
cpu_freq = (
    subprocess.check_output("grep 'cpu MHz' /proc/cpuinfo | head -1", shell=True)
    .decode()
    .split(":")[1]
    .strip()
)
print(f"The computation was performed on {cpu_name} operating at {cpu_freq}MHz.")


# In[9]:


import matplotlib.ticker as mticker

fig, ax = plt.subplots()
ax.set_title(
    rf"Surrogate evaluation time for $f_{{\rm{{start}}}} = {f_start} \rm{{Hz}}$, $f_{{\rm{{s}}}} = {1/delta_t:.0f} \rm{{Hz}}$"
)
ax.set_xlabel("Total mass ($M_\\odot$)")
ax.set_ylabel("Evaluation time (ms)")
ax.set_yscale("log")

alpha_scatter = 0.6
alpha_fill = 0.2
color = "dodgerblue"

for mass, eval_times in sur_eval_times_dict.items():
    eval_times = np.asarray(eval_times)

sur_eval_times_array = np.asarray(list(sur_eval_times_dict.values()))
sur_eval_times_median = np.median(sur_eval_times_array, axis=-1)
sur_speed_min = np.min(sur_eval_times_array, axis=-1)
sur_speed_max = np.max(sur_eval_times_array, axis=-1)

ax.plot(
    total_mass_array,
    sur_eval_times_median * 1.0e3,
    ls="--",
    color=color,
    marker=".",
    alpha=alpha_scatter,
    label="Median evaluation time",
)
ax.fill_between(
    total_mass_array,
    sur_speed_min * 1.0e3,
    sur_speed_max * 1.0e3,
    color=color,
    alpha=alpha_fill,
    label="Minimum/maximum evaluation time",
)


# To display only the first tick with a decimal entry
def mixed_formatter(x, pos):
    # pos = index of tick in visible ticks
    if pos == 0:
        return f"{x:.1f}"  # first tick: 1 decimal place
    else:
        return f"{x:.0f}"  # others: no decimals


ax.xaxis.set_major_formatter(mticker.FuncFormatter(mixed_formatter))
plt.xticks(total_mass_array)
plt.legend(frameon=False)
plt.tight_layout()

import subprocess

cpu_name = (
    subprocess.check_output("lscpu | grep 'Model name'", shell=True)
    .decode()
    .split(":")[1]
    .strip()
)
cpu_freq = (
    subprocess.check_output("grep 'cpu MHz' /proc/cpuinfo | head -1", shell=True)
    .decode()
    .split(":")[1]
    .strip()
)
print(f"The computation was performed on {cpu_name} operating at {cpu_freq}MHz.")


# In[16]:


import matplotlib.ticker as mticker

fig, ax = plt.subplots()
ax.set_title(
    rf"Surrogate evaluation time for $f_{{\rm{{start}}}} = {f_start} \rm{{Hz}}$, $f_{{\rm{{s}}}} = {1/delta_t:.0f} \rm{{Hz}}$"
)
ax.set_xlabel("Total mass ($M_\\odot$)")
ax.set_ylabel("Evaluation time (ms)")
ax.set_yscale("log")

alpha_scatter = 0.6
alpha_fill = 0.2
color = "dodgerblue"

for mass, eval_times in sur_eval_times_dict.items():
    eval_times = np.asarray(eval_times)

sur_eval_times_array = np.asarray(list(sur_eval_times_dict.values()))
sur_eval_times_median = np.median(sur_eval_times_array, axis=-1)
sur_speed_min = np.min(sur_eval_times_array, axis=-1)
sur_speed_max = np.max(sur_eval_times_array, axis=-1)

ax.plot(
    total_mass_array,
    sur_eval_times_median * 1.0e3,
    ls="--",
    color=color,
    marker=".",
    alpha=alpha_scatter,
    label="Median evaluation time",
)
ax.fill_between(
    total_mass_array,
    sur_speed_min * 1.0e3,
    sur_speed_max * 1.0e3,
    color=color,
    alpha=alpha_fill,
    label="Minimum/maximum evaluation time",
)


# To display only the first tick with a decimal entry
def mixed_formatter(x, pos):
    # pos = index of tick in visible ticks
    if pos == 0:
        return f"{x:.1f}"  # first tick: 1 decimal place
    else:
        return f"{x:.0f}"  # others: no decimals


ax.xaxis.set_major_formatter(mticker.FuncFormatter(mixed_formatter))
plt.xticks(total_mass_array)
plt.legend(frameon=False)
plt.tight_layout()

import subprocess

cpu_name = (
    subprocess.check_output("lscpu | grep 'Model name'", shell=True)
    .decode()
    .split(":")[1]
    .strip()
)
cpu_freq = (
    subprocess.check_output("grep 'cpu MHz' /proc/cpuinfo | head -1", shell=True)
    .decode()
    .split(":")[1]
    .strip()
)
print(f"The computation was performed on {cpu_name} operating at {cpu_freq}MHz.")


# <hr style="border: 2px solid #555; margin: 20px 0;">

# ## 4. The JAX backend
#
# Besides the default `NumPy` backend used so far, `InspiralESIGMASur` ships a [JAX](https://jax.readthedocs.io) backend (`esigmapy.surrogate.surrogate_jax`). It evaluates the *same* surrogate data through JAX-compiled kernels and reproduces the `NumPy` backend's waveforms to ~$10^{-11}$ relative accuracy, while unlocking three capabilities the `NumPy` backend cannot offer:
#
# 1. **Batched evaluations**: evaluate the waveform at *many* parameter-space points $(q, e_{\rm{ref}}, l_{\rm{ref}})$ in one vectorized call via `jax.vmap`, which amortizes the per-call overheads.
# 2. **Exact gradients**: differentiate the waveform (or any scalar built from it) with respect to the parameters via `jax.grad` — by automatic differentiation of the actual surrogate kernels, with no finite-difference stencils.
# 3. **Hardware portability**: the identical code runs on CPUs, GPUs and TPUs (more on this below).
#
# Two things to know before using it:
#
# - **Import order matters**: importing the backend enables JAX's 64-bit mode, which JAX fixes when the *first* JAX array of the process is created. Import `esigmapy.surrogate.surrogate_jax` before any other code that touches JAX.
# - **First calls compile**: the first evaluation of each JAX code path is traced and compiled by XLA (taking up to a few seconds); every later call reuses the compiled kernel. Never judge the JAX backend by its first call.
#
# The backend requires the JAX backend of `TPI` (`TPI_jax`); point the environment variable `TPI_JAX_PATH` to the `TPI` source tree if it is not importable already.

# ### 4.1 Constructing and evaluating the surrogate
#
# The JAX backend exposes the base surrogate classes directly. Construction loads the same data files as the `NumPy` backend and additionally moves them to the compute device.

# In[17]:


from esigmapy.surrogate.surrogate_jax import EccentricSurrogateJAX
from esigmapy.surrogate.surrogate import EccentricSurrogate
import jax
import jax.numpy as jnp

sur_jax = EccentricSurrogateJAX(
    ecc_data_dir=os.path.join(sur_data_dir, "ecc_sur_data"),
    circ_data_dir=os.path.join(sur_data_dir, "circ_sur_data"),
)
# The NumPy-backend surrogate object, used for comparisons below
sur_np = EccentricSurrogate(
    ecc_data_dir=os.path.join(sur_data_dir, "ecc_sur_data"),
    circ_data_dir=os.path.join(sur_data_dir, "circ_sur_data"),
)


# Single evaluations use the same call signature as the `NumPy` backend's base surrogate class: they return the time grid and the complex $(2,2)$ mode. Below we also overlay the `NumPy` backend result — the two are visually indistinguishable (we quantify the difference in section 4.6).

# In[18]:


M_tot = 10.0  # total mass (in solar masses)
q = 2.3
reference_eccentricity = 0.3
reference_mean_anomaly = 1.3
delta_t = 1 / 2**12
t_start = -2.0

t_jax, h_jax = sur_jax(
    M=M_tot,
    params=(q, reference_eccentricity, reference_mean_anomaly),
    delta_t=delta_t,
    t_start=t_start,
)
t_np, h_np = sur_np(
    M=M_tot,
    params=(q, reference_eccentricity, reference_mean_anomaly),
    delta_t=delta_t,
    t_start=t_start,
)

plt.figure(figsize=(10, 4))
plt.title(
    rf"$M={M_tot} M_\odot, q={q}, e_{{\rm{{ref}}}}={reference_eccentricity}, l_{{\rm{{ref}}}}={reference_mean_anomaly}$"
)
plt.plot(t_jax, h_jax.real, color="#D97706", label="JAX backend")
plt.plot(t_np, h_np.real, color="#1E90FF", ls="--", label="NumPy backend")
plt.xlabel(r"$t (s)$")
plt.ylabel(r"$\Re(h_{22})$")
plt.legend(frameon=False)
plt.tight_layout()


# ### 4.2 Batched evaluations over the parameter space
#
# The key entry point for the JAX-only features is `parameter_space_evaluator`. It fixes the time-grid configuration (`M`, `delta_t`, `t_start`/`t_end` or `times`) once on the host, and returns a *pure function* of the waveform parameters,
#
# ```python
# fn(q, e_ref, l_ref) -> (amp, phase)
# ```
#
# which composes with the standard JAX transformations: wrap it in `jax.jit` for fast repeated calls, in `jax.vmap` to evaluate whole parameter batches at once, and in `jax.grad` to differentiate. The complex mode is `amp * exp(-1j * phase)`.
#
# > **Note:** unlike the surrogate's `__call__`, `fn` performs *no* parameter-range checks (JAX-traced values cannot be inspected) and does not fall back to the circular surrogate at tiny $e_{\rm{ref}}$ — out-of-range inputs are silently extrapolated. Validate parameters with `sur_jax.check_param_range` first if in doubt.

# In[19]:


# Fix the time-grid configuration once...
t_grid, fn = sur_jax.parameter_space_evaluator(
    M=M_tot,
    delta_t=delta_t,
    t_start=t_start,
)

# ...and evaluate a whole batch of parameter-space points in ONE call.
q_batch = jnp.array([1.2, 2.3, 4.0, 5.5])
e_batch = jnp.array([0.05, 0.15, 0.3, 0.42])
l_batch = jnp.array([0.0, 1.3, 3.1, 5.0])

batched_fn = jax.jit(jax.vmap(fn))
amp_batch, phase_batch = batched_fn(q_batch, e_batch, l_batch)
print("batched amplitude array shape:", amp_batch.shape)  # (batch, num_samples)

fig, axs = plt.subplots(len(q_batch), sharex=True, figsize=(10, 8))
fig.suptitle(rf"One vmapped call: {len(q_batch)} waveforms at $M={M_tot}M_\odot$")
for i, ax in enumerate(axs):
    h_i = np.asarray(amp_batch[i]) * np.exp(-1j * np.asarray(phase_batch[i]))
    ax.plot(t_grid, h_i.real, color="#9467BD", lw=1.0)
    ax.set_ylabel(rf"$\Re(h_{{22}})$")
    ax.text(
        0.02,
        0.9,
        rf"$q={q_batch[i]:.1f}, e_{{\rm{{ref}}}}={e_batch[i]:.2f}, l_{{\rm{{ref}}}}={l_batch[i]:.1f}$",
        transform=ax.transAxes,
        va="top",
    )
plt.xlabel(r"$t (s)$")
plt.tight_layout()


# ### 4.3 Gradients in the parameter space
#
# Because `fn` is built entirely from differentiable JAX operations, `jax.grad` differentiates *through the surrogate* with respect to $(q, e_{\rm{ref}}, l_{\rm{ref}})$. This is exact automatic differentiation of the surrogate's own kernels (splines, empirical-interpolant contractions, interpolations) — no finite-difference step-size tuning, no extra surrogate evaluations.
#
# As an example we differentiate a simple scalar functional of the waveform, its squared amplitude norm $\mathcal{E} = \sum_i A_i^2$, and check the eccentricity derivative against a central finite difference. Gradients like these are what gradient-based samplers (Hamiltonian Monte Carlo, variational inference) and Fisher-matrix codes need.

# In[20]:


def energy(q, e_ref, l_ref):
    amp, _ = fn(q, e_ref, l_ref)
    return jnp.sum(amp**2)


# d(energy)/d(e_ref), exact via automatic differentiation
denergy_deref = jax.jit(jax.grad(energy, argnums=1))
g = float(denergy_deref(2.3, 0.3, 1.3))

# Finite-difference cross-check. NOTE: the functional is only piecewise smooth
# (the surrogate interpolates tables linearly in time), so the FD step must be
# small enough to stay within one smooth piece.
eps = 1e-7
fd = (float(energy(2.3, 0.3 + eps, 1.3)) - float(energy(2.3, 0.3 - eps, 1.3))) / (
    2 * eps
)
print(f"jax.grad : {g:.10e}")
print(f"finite difference: {fd:.10e}   (relative difference {abs(g - fd) / abs(fd):.1e})")

# Transformations compose: the gradient at a whole batch of points, in one call
batch_grads = jax.vmap(jax.grad(energy, argnums=(0, 1, 2)))(q_batch, e_batch, l_batch)
print("batched (dE/dq, dE/de_ref, dE/dl_ref) shapes:", [b.shape for b in batch_grads])


# ### 4.4 Gradients along time
#
# The waveform is returned as sampled arrays, so time derivatives are also available by differentiating through an interpolant of the output — e.g. the instantaneous GW frequency $f_{\rm{GW}}(t) = \frac{1}{2\pi} \frac{\mathrm{d}\phi}{\mathrm{d}t}$ from the phase. Here we differentiate `jnp.interp` over the sampled phase with `jax.grad`, vectorized over the whole grid with `jax.vmap`.
#
# > **Note:** the sampled series is piecewise linear in time, so this derivative is piecewise constant between samples; at the surrogate's dense output sampling this is indistinguishable from a smooth curve on any plot, but keep it in mind if you need higher-order time derivatives.

# In[21]:


amp_1, phase_1 = jax.jit(fn)(2.3, 0.3, 1.3)
t_grid_j = jnp.asarray(t_grid)


def phase_at(t):
    return jnp.interp(t, t_grid_j, phase_1)


# dphi/dt at every grid point (evaluated slightly inside the ends), in one call
t_eval = t_grid_j[1:-1]
f_gw = jax.vmap(jax.grad(phase_at))(t_eval) / (2 * np.pi)

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 6))
ax1.plot(t_grid, np.asarray(amp_1), color="#D97706")
ax1.set_ylabel(r"$A_{22}$")
ax2.plot(np.asarray(t_eval), np.asarray(f_gw), color="#D97706")
ax2.set_ylabel(r"$f_{\rm{GW}} \, (\rm{Hz})$")
plt.xlabel(r"$t (s)$")
fig.suptitle(
    r"Amplitude and instantaneous GW frequency via $\frac{1}{2\pi}\mathrm{d}\phi/\mathrm{d}t$ (note the eccentric oscillations)"
)
plt.tight_layout()


# ### 4.5 Evaluation speed: single, batched, and vs the NumPy backend
#
# We now repeat the timing experiment of section 3.4 with three configurations:
#
# 1. **JAX, single evaluation**: one jitted `fn(q, e_ref, l_ref)` call at a time;
# 2. **JAX, batched**: `jax.vmap(fn)` over all parameter draws of a mass in one call, reported *per waveform*;
# 3. **NumPy backend, single evaluation** — the baseline from section 3.4.
#
# As before, waveforms start at a GW frequency of $15 \rm{Hz}$ and are sampled at $4096 \rm{Hz}$. Compilation (warmup) is excluded from all JAX timings: each jitted function is called once before the timed repetitions, and every timing waits for the device with `jax.block_until_ready`. The batched configuration is where JAX shines — the fixed Python-to-XLA dispatch cost and the kernel launch overheads are paid once per *batch* instead of once per *waveform*, and this is exactly the shape of the workload in parameter estimation.
#
# > **Note:** the following cell may take a few minutes to run.

# In[22]:


import time
import statistics
from esigmapy.utils import x_from_f22

f_start = 15.0
delta_t = 1 / 2**12
num_evals_per_mass = 25
NLOOP = 10

# Same mass grid as section 3.4 (restricted to masses whose full surrogate
# length reaches 15 Hz; see section 3.3 for M_min_at_fcompute)
bench_mass_array = np.asarray([20, 30, 40, 50, 60, 80, 100])
bench_mass_array = bench_mass_array[bench_mass_array > M_min_at_fcompute]
bench_mass_array = np.r_[M_min_at_fcompute, bench_mass_array]


def median_time(callable_fn, nloop=NLOOP):
    """Median wall-clock seconds per call; the caller is responsible for
    warming the jit cache first. block_until_ready is applied by callable_fn."""
    samples = []
    for _ in range(nloop):
        tic = time.perf_counter()
        callable_fn()
        samples.append(time.perf_counter() - tic)
    return statistics.median(samples)


def t_start_for_mass(M):
    """Starting time at which a representative waveform reaches f_start."""
    t, orb_vars, _ = esigmasur.get_inspiral_esigmasur_modes(
        mass1=3.0 * M / 4.0,
        mass2=M / 4.0,
        reference_eccentricity=0.3,
        reference_mean_anomaly=np.pi,
        delta_t=delta_t,
        t_start=None,
        return_pycbc_timeseries=False,
        return_orbital_params=["x"],
    )
    idx = np.argmax(orb_vars["x"] >= x_from_f22(f_start, M=M))
    return t[idx]


rng = np.random.default_rng(seed=37)
jax_single_times = []
jax_batched_times = []
numpy_single_times = []

for M in bench_mass_array:
    q_arr = rng.uniform(1, 6, size=num_evals_per_mass)
    # stay above the tiny-e_ref circular-fallback region, which
    # parameter_space_evaluator does not handle (see the note in 4.2)
    e_arr = rng.uniform(0.01, 0.43, size=num_evals_per_mass)
    l_arr = rng.uniform(0, 2 * np.pi, size=num_evals_per_mass)
    ts = t_start_for_mass(M)

    _, fn_M = sur_jax.parameter_space_evaluator(M=M, delta_t=delta_t, t_start=ts)
    fn_single = jax.jit(fn_M)
    fn_batch = jax.jit(jax.vmap(fn_M))

    # -- JAX single (median per-call over the draws; jit warmed by the first draw)
    jax.block_until_ready(fn_single(q_arr[0], e_arr[0], l_arr[0]))
    per_call = [
        median_time(
            lambda i=i: jax.block_until_ready(fn_single(q_arr[i], e_arr[i], l_arr[i])),
            nloop=3,
        )
        for i in range(num_evals_per_mass)
    ]
    jax_single_times.append(np.median(per_call))

    # -- JAX batched (one vmapped call for all draws, reported per waveform)
    qj, ej_, lj = jnp.asarray(q_arr), jnp.asarray(e_arr), jnp.asarray(l_arr)
    jax.block_until_ready(fn_batch(qj, ej_, lj))
    batch_time = median_time(lambda: jax.block_until_ready(fn_batch(qj, ej_, lj)))
    jax_batched_times.append(batch_time / num_evals_per_mass)

    # -- NumPy backend single
    per_call = [
        median_time(
            lambda i=i: sur_np(
                M=M, params=(q_arr[i], e_arr[i], l_arr[i]), delta_t=delta_t, t_start=ts
            ),
            nloop=3,
        )
        for i in range(num_evals_per_mass)
    ]
    numpy_single_times.append(np.median(per_call))

    print(f"Timing completed for M = {M:.2f} M_Sun")


# In[23]:


import matplotlib.ticker as mticker

fig, ax = plt.subplots()
ax.set_title(
    rf"Waveform evaluation cost for $f_{{\rm{{start}}}} = {f_start} \rm{{Hz}}$, $f_{{\rm{{s}}}} = {1/delta_t:.0f} \rm{{Hz}}$"
)
ax.set_xlabel("Total mass ($M_\\odot$)")
ax.set_ylabel("Evaluation time per waveform (ms)")
ax.set_yscale("log")

ax.plot(
    bench_mass_array,
    np.asarray(numpy_single_times) * 1.0e3,
    color="#1E90FF",
    ls="--",
    marker="o",
    label="NumPy backend, single evaluation",
)
ax.plot(
    bench_mass_array,
    np.asarray(jax_single_times) * 1.0e3,
    color="#D97706",
    ls="-",
    marker="s",
    label="JAX backend, single evaluation",
)
ax.plot(
    bench_mass_array,
    np.asarray(jax_batched_times) * 1.0e3,
    color="#9467BD",
    ls="-",
    marker="^",
    label=rf"JAX backend, batched ($\times${num_evals_per_mass}), per waveform",
)

ax.xaxis.set_major_formatter(mticker.FuncFormatter(mixed_formatter))
plt.xticks(bench_mass_array)
plt.legend(frameon=False)
plt.tight_layout()

cpu_name = (
    subprocess.check_output("lscpu | grep 'Model name'", shell=True)
    .decode()
    .split(":")[1]
    .strip()
)
print(f"The computation was performed on {cpu_name} (CPU).")


# ### A note on hardware: CPUs, GPUs and TPUs
#
# Everything in this section ran on a CPU, but none of the code above is CPU-specific: JAX compiles the same Python to whatever accelerator backend is installed. With a CUDA-enabled `jaxlib` (`pip install "jax[cuda12]"`) and an NVIDIA GPU visible to the process, the surrogate data is placed on the GPU at construction and all evaluations, batches and gradients execute there — no changes to user code.
#
# GPUs pay off most for exactly the workloads the JAX backend targets: **large batched evaluations** (and batched gradients), where thousands of waveforms are computed by one massively parallel kernel launch, can achieve substantially better per-waveform costs than the CPU curves above. Single-waveform latency, by contrast, is typically *not* better on a GPU — the fixed host-to-device dispatch overhead dominates a single short kernel.
#
# Two GPU-specific caveats:
#
# - The surrogate requires **64-bit floats** (`TPI_jax` enables JAX's x64 mode on import). Consumer/gaming GPUs execute float64 at a small fraction of their float32 rate; data-center GPUs (A100/H100 class) have full-rate float64 and are the appropriate targets for production batches.
# - Very large batches of long waveforms are memory-bound: the batch of output arrays must fit in device memory; chunk the batch (e.g. with `jax.lax.map` over vmapped chunks) if needed.

# ### 4.6 Cost of parameter-space gradients
#
# The same batching logic applies to gradients. We time the gradient of the scalar functional from section 4.3, $\nabla_{(q, e_{\rm{ref}}, l_{\rm{ref}})} \sum_i A_i^2$, as (1) a single jitted `jax.grad` evaluation and (2) a `jax.vmap`-batched gradient over all draws of a mass, per point. A reverse-mode gradient evaluation costs a small constant multiple of the forward evaluation — compare with the previous plot — and batching amortizes its overheads the same way.
#
# > **Note:** the following cell may take a few minutes to run.

# In[24]:


jax_grad_single_times = []
jax_grad_batched_times = []

rng = np.random.default_rng(seed=37)
for M in bench_mass_array:
    q_arr = rng.uniform(1, 6, size=num_evals_per_mass)
    e_arr = rng.uniform(0.01, 0.43, size=num_evals_per_mass)
    l_arr = rng.uniform(0, 2 * np.pi, size=num_evals_per_mass)
    ts = t_start_for_mass(M)

    _, fn_M = sur_jax.parameter_space_evaluator(M=M, delta_t=delta_t, t_start=ts)

    def energy_M(q, e_ref, l_ref, fn_M=fn_M):
        amp, _ = fn_M(q, e_ref, l_ref)
        return jnp.sum(amp**2)

    grad_single = jax.jit(jax.grad(energy_M, argnums=(0, 1, 2)))
    grad_batch = jax.jit(jax.vmap(jax.grad(energy_M, argnums=(0, 1, 2))))

    jax.block_until_ready(grad_single(q_arr[0], e_arr[0], l_arr[0]))
    per_call = [
        median_time(
            lambda i=i: jax.block_until_ready(grad_single(q_arr[i], e_arr[i], l_arr[i])),
            nloop=3,
        )
        for i in range(num_evals_per_mass)
    ]
    jax_grad_single_times.append(np.median(per_call))

    qj, ej_, lj = jnp.asarray(q_arr), jnp.asarray(e_arr), jnp.asarray(l_arr)
    jax.block_until_ready(grad_batch(qj, ej_, lj))
    batch_time = median_time(lambda: jax.block_until_ready(grad_batch(qj, ej_, lj)))
    jax_grad_batched_times.append(batch_time / num_evals_per_mass)

    print(f"Gradient timing completed for M = {M:.2f} M_Sun")


# In[25]:


fig, ax = plt.subplots()
ax.set_title(
    rf"Parameter-space gradient cost for $f_{{\rm{{start}}}} = {f_start} \rm{{Hz}}$, $f_{{\rm{{s}}}} = {1/delta_t:.0f} \rm{{Hz}}$"
)
ax.set_xlabel("Total mass ($M_\\odot$)")
ax.set_ylabel("Gradient evaluation time per point (ms)")
ax.set_yscale("log")

ax.plot(
    bench_mass_array,
    np.asarray(jax_grad_single_times) * 1.0e3,
    color="#D97706",
    ls="-",
    marker="s",
    label=r"$\nabla_{(q, e_{\rm{ref}}, l_{\rm{ref}})}$, single evaluation",
)
ax.plot(
    bench_mass_array,
    np.asarray(jax_grad_batched_times) * 1.0e3,
    color="#9467BD",
    ls="-",
    marker="^",
    label=rf"$\nabla_{{(q, e_{{\rm{{ref}}}}, l_{{\rm{{ref}}}})}}$, batched ($\times${num_evals_per_mass}), per point",
)

ax.xaxis.set_major_formatter(mticker.FuncFormatter(mixed_formatter))
plt.xticks(bench_mass_array)
plt.legend(frameon=False)
plt.tight_layout()


# ### 4.7 Are the two backends equal? Mismatch distribution
#
# Finally, we quantify the agreement between the two backends with the standard waveform-comparison metric: the **mismatch** $1 - \mathcal{O}$, where $\mathcal{O}$ is the normalized (time-domain, flat-noise) overlap
#
# $$ \mathcal{O}(h_a, h_b) = \frac{|\langle h_a, h_b \rangle|}{\sqrt{\langle h_a, h_a \rangle \langle h_b, h_b \rangle}}. $$
#
# We draw random parameter-space points, generate the waveform with both backends on the identical time grid, and histogram the mismatches. The values sit around $10^{-22}\text{–}10^{-16}$ — at (or below) double-precision resolution of the overlap, i.e. the two backends are equal for every practical purpose (waveform differences are $\sim 10^{-11}$ in relative amplitude, and the mismatch is quadratic in the difference).
#
# > **Note:** the following cell may take a few minutes to run.

# In[26]:


num_mismatch_draws = 200
M_mm = 10.0
t_start_mm = -10.0

rng = np.random.default_rng(seed=42)
q_arr = rng.uniform(1, 6, size=num_mismatch_draws)
e_arr = rng.uniform(0.01, 0.43, size=num_mismatch_draws)
l_arr = rng.uniform(0, 2 * np.pi, size=num_mismatch_draws)


def mismatch(h_a, h_b):
    overlap = np.abs(np.vdot(h_a, h_b)) / np.sqrt(
        np.vdot(h_a, h_a).real * np.vdot(h_b, h_b).real
    )
    return 1.0 - overlap


_, fn_mm = sur_jax.parameter_space_evaluator(M=M_mm, delta_t=delta_t, t_start=t_start_mm)
fn_mm = jax.jit(fn_mm)

mismatches = []
for i in range(num_mismatch_draws):
    amp_j, phase_j = fn_mm(q_arr[i], e_arr[i], l_arr[i])
    h_j = np.asarray(amp_j) * np.exp(-1j * np.asarray(phase_j))
    _, h_n = sur_np(
        M=M_mm,
        params=(q_arr[i], e_arr[i], l_arr[i]),
        delta_t=delta_t,
        t_start=t_start_mm,
    )
    mismatches.append(mismatch(h_j, h_n))

mismatches = np.asarray(mismatches)
# mismatches can round to exactly 0 at double precision; floor them at 1e-22
# so they remain visible on the logarithmic axis
mismatches = np.maximum(mismatches, 1e-22)
print(f"median mismatch: {np.median(mismatches):.2e}, max: {mismatches.max():.2e}")


# In[27]:


fig, ax = plt.subplots()
ax.set_title(
    rf"JAX vs NumPy backend mismatches: {num_mismatch_draws} random draws at $M={M_mm:.0f}M_\odot$"
)
bins = np.logspace(-22, -12, 41)
ax.hist(mismatches, bins=bins, color="#1E90FF", edgecolor="white", linewidth=0.5)
ax.set_xscale("log")
ax.set_xlabel(r"mismatch $1-\mathcal{O}$ (values at $10^{-22}$ are floored)")
ax.set_ylabel("Number of draws")
ax.axvline(1e-15, color="0.4", ls=":", lw=1)
ax.text(
    1.2e-15,
    ax.get_ylim()[1] * 0.9 if ax.get_ylim()[1] > 0 else 1,
    "double-precision scale",
    rotation=90,
    va="top",
    color="0.4",
    fontsize=10,
)
plt.tight_layout()


# In[ ]:




