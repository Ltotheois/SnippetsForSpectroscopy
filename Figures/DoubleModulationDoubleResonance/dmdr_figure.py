#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Luis Bonah
# Description: Plots Detector Signal and Demodulation of DM-DR


##
## Modules
##

import sys
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy import special


##
## Settings
##

kwargs_figure = {
	"figsize" :						(16/2.54, 10/2.54),
	"dpi" :							300,
}
kwargs_plot = {
	"linewidth" :					1,
}

kwargs_font = {
	'size' :						6,
}
matplotlib.rc('font', **kwargs_font)

settings = {
	"splitting" 		:	3.0,
	"width_gauss"		:	0.5,
	"width_lorentz"		:	0.5,
	"amp"				:	np.float64(1/(2*0.20870928052036772)),
	"fm_freq"			:	np.float64(0.10),
	"am_freq"			:	np.float64(0.02),
	"fm_phase"			:	np.pi/2,
	"am_phase"			:	0,
	"hub"				:	1.5,
	"cp"				:	0.0,
	"duration"			:	100,
	"frames"			:	1000,
	"inline_sig_start"	:	4,
	"inline_fm_start"	:	1.2,
	"pixels_spectrum"	:	1000,
	"spectrum_width"	:	20,
	"spectrum_height"	:	2,
}

arrow_dict = {
	"head_width" :		0.1,
	"head_length" :		0.1,
	"length_includes_head" : True,
	"width" :			0.002,
	"zorder" :			100,
	"color" :			"#737373",
}

fd = {
	"color" : "#737373",
	"size"	: 4,
}

color_sin  = "#0336FF"
color_off  = "#6ebeff"
color_on   = "#FF0266"
color_help = "#c2c2c2"
color_0    = "#878787"
color_pos  = "#CCCCCC"
color_neg  = "#89919A"
color_dmdr = "#000000"

##
## Helper Function
##

def shared_labels(fig, gsp, xlabel, ylabel):
	ax = fig.add_subplot(gsp, frameon=False)
	for label in ("top", "bottom", "right", "left"):
		ax.spines[label].set_color('none')
	ax.tick_params(axis="both", labelcolor='#fff0', top=False, bottom=False, left=False, right=False, zorder=-100)
	ax.set_xlabel(xlabel, labelpad = 0)
	ax.set_ylabel(ylabel, labelpad = 5)
	return(ax)

def lineshape_off(xs):
	ys = special.voigt_profile(xs, settings["width_gauss"], settings["width_lorentz"])
	return(ys*settings["amp"])
	
def lineshape_on(xs):
	ys = 0.5*(special.voigt_profile(xs-settings["splitting"]/2, settings["width_gauss"], settings["width_lorentz"])+special.voigt_profile(xs+settings["splitting"]/2, settings["width_gauss"], settings["width_lorentz"]))
	return(ys*settings["amp"])

def integrate_parts(signal, parts):
	tmp_integral = []
	tmp_output = []
	tmp_current_i = 0
	for i, y in zip(parts, signal):
		if i != tmp_current_i:
			mean = sum(tmp_integral)/len(tmp_integral)
			tmp_output.extend([mean]*len(tmp_integral))
			tmp_integral = []
			tmp_current_i = i
		tmp_integral.append(y)
	tmp_output.extend([np.nan]*len(tmp_integral))
	return(np.array(tmp_output))

def make_collection(ax, xs, ys, colors, kwargs_plot, ind=None):
	if ind == None:
		ind = len(xs)-1
	if ind <= 0:
		segs = (((0,0),(0, 0)) for i in range(0))
	else:
		segs = (((xs[i], ys[i]),(xs[i+1], ys[i+1])) for i in range(ind))
	coll = matplotlib.collections.LineCollection(segs, color=colors, **kwargs_plot)
	ref = ax.add_collection(coll)
	return(ref)

##
## Basic Setup
##

labels = [1, "h1", 2, 3, 4, 6, 7, "h2", 8]
ratios = [1 if str(labels[i]).isnumeric() else 0.85 for i in range(len(labels))]
ratios[3] = 0.5
ratios[6] = 0.5

ncols = 11
gs = matplotlib.gridspec.GridSpec(len(labels), ncols, height_ratios = ratios, hspace=0)

fig = plt.figure(**kwargs_figure)
axs = {}

for x, label in enumerate(labels[1:]):
	if label == 8:
		axs[label] = fig.add_subplot(gs[x+1,3:8])
	else:
		axs[label] = fig.add_subplot(gs[x+1,:])
axs["l"] = fig.add_subplot(gs[0,:int((ncols-1)/2)])
axs["r"] = fig.add_subplot(gs[0,int((ncols+1)/2):])

for x in [1,2,3,4,5]:
	axs[labels[x]].get_xaxis().set_visible(False)


##
## Prepare Data
##

spectrum_xs = np.linspace(-settings["spectrum_width"]/2, settings["spectrum_width"]/2, settings["pixels_spectrum"])
spectrum_ys_off = lineshape_off(spectrum_xs)
spectrum_ys_on  = lineshape_on(spectrum_xs)

time_xs = np.linspace(0, settings["duration"], settings["frames"]+1)
freq_xs = settings["cp"] + settings["hub"]*np.sin(2*np.pi*time_xs*settings["fm_freq"])
am_sine_ys = np.sin(2*np.pi*settings["am_freq"]*time_xs + settings["am_phase"])
fm2_sine_ys = np.sin(2*2*np.pi*settings["fm_freq"]*time_xs + settings["fm_phase"])

detector_ys = np.where(am_sine_ys >= 0, lineshape_off(freq_xs), lineshape_on(freq_xs))
mixed_ys = fm2_sine_ys * detector_ys

int_ind = time_xs // (1/settings["am_freq"]/2)
lockin_ys = integrate_parts(mixed_ys, int_ind)
mixed2_ys = mixed_ys*am_sine_ys

colors = np.where(am_sine_ys >= 0, color_off, color_on)
max_signal = max(np.nanmax(spectrum_ys_off), np.nanmax(spectrum_ys_on))*1.2

# Inline FM Sine
tmp_min = max_signal
tmp_max = settings["spectrum_height"]
inline_sine_xs = tmp_min+time_xs*(tmp_max-tmp_min)/settings["duration"]

# Inline Detector Signal
tmp_min = settings["inline_sig_start"]
tmp_max = settings["spectrum_width"]/2
inline_det_xs = tmp_min+time_xs*(tmp_max-tmp_min)/settings["duration"]

# Resulting Point in DMDR Spectrum
resulting_point = (settings["cp"], np.nanmean(mixed2_ys))


##
## Uncomment to test initial data
##

# axs[1].plot(spectrum_xs, spectrum_ys_on)
# axs[1].plot(spectrum_xs, spectrum_ys_off)
# axs[2].plot(time_xs, detector_ys)
# axs[3].plot(time_xs, fm2_sine_ys)
# axs[4].plot(time_xs, mixed_ys)
# axs[5].plot(time_xs, lockin_ys)
# axs[6].plot(time_xs, am_sine_ys)
# axs[7].plot(time_xs, mixed2_ys)
# plt.show()

##
## Draw Initial Plots
##

_ = shared_labels(fig, gs[2:7, :], "Time /a.u.", "Int. /a.u.")


# --- First Row ---
for ax, color, func in zip([axs["l"], axs["r"]], [color_off, color_on], [lineshape_off, lineshape_on]):
	ax.plot(freq_xs, inline_sine_xs, color = "#c2c2c2"+"82", linewidth = 0.5)
	ax.plot(inline_det_xs, func(freq_xs), color = color+"82", linewidth = 0.5)
	ax.plot(spectrum_xs, func(spectrum_xs), color = color)
	
	lxs = [0, 0.5*settings["hub"], settings["hub"]]
	lys = [func(x) for x in lxs]
	ax.vlines(lxs, lys, +100, color = "#c2c2c2"+"82", linewidth = 0.5)
	ax.hlines(lys, lxs, +100, color = "#c2c2c2"+"82", linewidth = 0.5)
	ax.arrow(0, 1.2, 0, 0.735, **arrow_dict)
	ax.arrow(4, (max(lys)+min(lys))/2, 5.9, 0, **arrow_dict)
	ax.text(0, 1.95, "  $t$", ha="left", va="top", transform=ax.transData, **fd)
	ax.text(9.9, (max(lys)+min(lys))/2+0.05, "$t$ ", ha="right", va="bottom", transform=ax.transData, **fd)
	ax.text(+settings["hub"], 1.95, " FM ", ha="left", va="top", transform=ax.transData, size=4, color="#c2c2c2")
	ax.text(+4, max_signal*1.2, "Signal at Detector", ha="left", va="center", transform=ax.transData, size=4, color=color+"82")

	ax.set_xlim(-settings["spectrum_width"]/2, settings["spectrum_width"]/2)
	ax.set_ylim(0, settings["spectrum_height"])
	ax.set_xlabel("Probe Frequency /a.u.")
	
	ax.set_xticks(np.linspace(-settings["spectrum_width"]/2, settings["spectrum_width"]/2, 5))

axs["l"].get_shared_y_axes().join(axs["l"], axs["r"])
axs["l"].get_shared_x_axes().join(axs["l"], axs["r"])

axs["l"].set_ylabel("Int. /a.u.", labelpad = 10)
axs["r"].text(1.0, 0.5, "  1) Spectrum", ha="left", va="center", transform=axs["r"].transAxes)


# --- Hidden Row ---

axs["h1"].axis("off")


# --- Second Row ---

# Detector Signal
plot_detector = make_collection(axs[2], time_xs, detector_ys, colors, kwargs_plot)

# Styling 
axs[2].set_xlim(0, settings["duration"])
axs[2].set_ylim(0, max_signal)
axs[2].text(1.0, 0.5, "  2) Signal at Detector $D(t)$", ha="left", va="center", transform=axs[2].transAxes)
axs[2].axhline(y=0, color=color_0, linewidth = 0.25)


# --- Third Row ---

# 2f Signal
axs[3].plot(time_xs, fm2_sine_ys, color=color_sin, **kwargs_plot)

# Styling
axs[3].get_yaxis().set_visible(False)
axs[3].set_xlim(0, settings["duration"])
axs[3].set_ylim((-1.8, 1.8))
axs[3].axhline(y=0, color=color_0, linewidth = 0.25)
axs[3].text(1.0, 0.5, r"  3) $f(t)=\sin(2\cdot 2\pi\cdot f\cdot t + \phi)$", ha="left", va="center", transform=axs[3].transAxes)


# --- Fourth Row ---

# Mixed Signal
plot_mixed = make_collection(axs[4], time_xs, mixed_ys, colors, kwargs_plot)

# # Filling
# fill1 = axs[4].fill_between(time_xs, mixed_ys, mixed_ys*0, mixed_ys>0, color=color_pos)
# fill2 = axs[4].fill_between(time_xs, mixed_ys*0, mixed_ys, mixed_ys<0, color=color_neg)

# Styling
axs[4].set_xlim(0, settings["duration"])
axs[4].set_ylim(np.nanmin(mixed_ys)*1.2 if np.nanmin(mixed_ys) < 0 else -1, np.nanmax(mixed_ys)*1.2 if np.nanmax(mixed_ys) > 0 else 1)
axs[4].axhline(y=0, color='#878787', linewidth = 0.25)
axs[4].text(1.0, 0.5, r"  4) $2f$ demod.: $S(t)=f(t)\cdot D(t)$", ha="left", va="center", transform=axs[4].transAxes)


# --- Sixth Row ---

# AM Signal
axs[6].plot(time_xs, am_sine_ys, color=color_sin, **kwargs_plot)

# Styling
axs[6].get_yaxis().set_visible(False)
axs[6].set_xlim(0, settings["duration"])
axs[6].set_ylim((-1.8, 1.8))
axs[6].axhline(y=0, color=color_0, linewidth = 0.25)
axs[6].text(1.0, 0.5, r"  5) $f'(t)=\sin(2\pi\cdot f'\cdot t + \varphi)$", ha="left", va="center", transform=axs[6].transAxes)

for j in range(0, int(settings["duration"]*settings["am_freq"])):
	va = "bottom"
	text = "+"
	axs[6].text((j+0.25)/settings["am_freq"], 0, text, ha="center", va=va, transform=axs[6].transData, color="#737373", size=9)
	axs[6].axvline((j)/settings["am_freq"], color="#737373", linewidth=0.2)
	
	va = "top"
	text = "–"
	axs[6].text((j+0.75)/settings["am_freq"], 0, text, ha="center", va=va, transform=axs[6].transData, color="#737373", size=9)
	axs[6].axvline((j+0.5)/settings["am_freq"], color="#737373", linewidth=0.2)

# --- Seventh Row ---

# Mixed 2 Lockin Signal
plot_mixed2 = make_collection(axs[7], time_xs, mixed2_ys, colors, kwargs_plot)

# Filling
fill3 = axs[7].fill_between(time_xs, mixed2_ys, mixed2_ys*0, mixed2_ys>0, color=color_pos)
fill4 = axs[7].fill_between(time_xs, mixed2_ys*0, mixed2_ys, mixed2_ys<0, color=color_neg)

# Styling
axs[7].set_xlim(0, settings["duration"])
axs[7].set_ylim(np.nanmin(mixed2_ys)*1.2 if np.nanmin(mixed2_ys) < 0 else 0, np.nanmax(mixed2_ys)*1.2 if np.nanmax(mixed2_ys) > 0 else 1)
axs[7].axhline(y=0, color=color_0, linewidth = 0.25)
axs[7].text(1.0, 0.5, r"  6) 1$f'$ demod.: $f'(t)\cdot S(t)$", ha="left", va="center", transform=axs[7].transAxes)

axs[7].set_xticks(np.linspace(0, settings["duration"], 5))
axs[7].set_xticklabels(f"{x*settings['am_freq']}" for x in np.linspace(0, settings["duration"], 5))

# --- Hidden Row ---

axs["h2"].axis("off")


# --- Eigth Row ---

# DMDR Signal
tmp_ys = np.gradient(np.gradient(spectrum_ys_on))-np.gradient(np.gradient(spectrum_ys_off))
tmp_ys = tmp_ys/tmp_ys[int(len(tmp_ys)/2)]*np.nanmean(mixed2_ys)
axs[8].plot(spectrum_xs, tmp_ys, color=color_dmdr+"40", **kwargs_plot)
plot_dmdr = axs[8].scatter(*resulting_point, color=color_dmdr, s=12)

# Styling
axs[8].set_xlim(-settings["spectrum_width"]/2, settings["spectrum_width"]/2)
axs[8].set_ylim(np.nanmin(tmp_ys)*1.2, np.nanmax(tmp_ys)*1.2)
axs[8].set_xlabel("Probe Frequency /a.u.")
axs[8].set_ylabel("Int. /a.u.", labelpad = 10)
axs[8].text(1.0, 0.5, "  7) DM-DR Spectrum", ha="left", va="center", transform=axs[8].transAxes)

axs[8].set_xticks(np.linspace(-settings["spectrum_width"]/2, settings["spectrum_width"]/2, 5))
yt = (0, 0.1)
axs[8].set_yticks(yt)
axs[8].set_yticklabels([f"{x*10:.0f}" for x in yt])


# --- Connections ---
conn_kwargs = {
	"axesA"		: axs[7],
	"axesB"		: axs[8],
	"coordsA"	: "data",
	"coordsB"	: "data",
	"arrowstyle": "-",
	"color"		: "#878787",
}
	
patch = matplotlib.patches.ConnectionPatch((0, 0), resulting_point, **conn_kwargs)
axs[8].add_artist(patch)

patch = matplotlib.patches.ConnectionPatch((settings["duration"], 0), resulting_point, **conn_kwargs)
axs[8].add_artist(patch)


##
## Main Work
##
fig.subplots_adjust(right = 0.8, top = 0.96)

plt.savefig("dmdr_figure.pdf")
# plt.show()
plt.close()