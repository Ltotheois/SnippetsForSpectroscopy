#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Luis Bonah, Oliver Zingsheim, Milan Roska
# Description: Plots different FM Amplitudes and Resulting Lineshapes


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import special

##
## Setup
##
kwargs_figure = {
	"figsize" :						(14/2.54, 10/2.54),
	"dpi" :							300,
}
kwargs_font = {
	'size' :						10,
}
kwargs_plot = {
	"linewidth" :					1,
}

matplotlib.rc('font', **kwargs_font)

settings = {
	"width_gauss"		:	0.5,
	"width_lorentz"		:	0.5,
	"amp"				:	np.float64(1/(2*0.20870928052036772)),
	"fm_freq"			:	np.float64(0.06),
	"fm1_phase"			:	0,
	"fm2_phase"			:	np.pi/2,
	"hub"				:	1.5,
	"cp"				:	1.5,
	"points_time"		:	100,
	"frames"			:	1000,
	"inline_sig_start"	:	4,
	"points_spectrum"	:	1000,
	"spectrum_width"	:	20,
	"spectrum_height"	:	3,
}


color_fm   = "#b0b0b0"
color_result = "#6ebeff"
color_result_2f = "#0336FF"
color_spectrum = "#FF0266"

fd = {
	"color" : "#666666",
	"size"	: 5,
}
arrow_dict = {
	"head_width" :		0.1,
	"head_length" :		0.1,
	"length_includes_head" : True,
	"width" :			0.001,
	"zorder" :			100,
	"color" :			"#666666",
	"linewidth":		0.5
}

##
## Functions
##
def shared_labels(fig, gsp, xlabel, ylabel):
	ax = fig.add_subplot(gsp, frameon=False)
	for label in ("top", "bottom", "right", "left"):
		ax.spines[label].set_color('none')
	ax.tick_params(axis="both", labelcolor='#fff0', top=False, bottom=False, left=False, right=False, zorder=-100)
	ax.set_xlabel(xlabel, labelpad = 5)
	ax.set_ylabel(ylabel, labelpad = 5)
	return(ax)

def lineshape(xs):	
	ys = special.voigt_profile(xs, settings["width_gauss"], settings["width_lorentz"])
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

def spectrum_dm_modulation(offsets, lineshape, freq_xs, sine):
	out = []
	for offset in offsets:
		tmp = lineshape(freq_xs+offset-settings["cp"])*sine
		res = np.sum(tmp)/len(tmp)
		out.append(res)
	return(np.array(out))


##
## Gridspecs
##
gs = GridSpec(3, 3, hspace=0, wspace=0)
fig = plt.figure(**kwargs_figure)

axs = {}
axs[0] = fig.add_subplot(gs[0, 0])
axs[10] = fig.add_subplot(gs[0, 1])
axs[20] = fig.add_subplot(gs[0, 2])
axs[1] = fig.add_subplot(gs[1, 0])
axs[11] = fig.add_subplot(gs[1, 1])
axs[21] = fig.add_subplot(gs[1, 2])
axs[2] = fig.add_subplot(gs[2, 0])
axs[12] = fig.add_subplot(gs[2, 1])
axs[22] = fig.add_subplot(gs[2, 2])

for x in [0, 1, 10, 11, 20, 21]:
	axs[x].get_xaxis().set_visible(False)


##
## Filling Plots
##
hubs = (0.5, 1.5, 3)
# hubs = (0.5, 1.89, 3)
res_range_1f = [0, 0]
res_range_2f = [0, 0]

for i, hub in enumerate(hubs):
	settings["hub"] = hub
	
	# Basic Arrays
	time_xs = np.linspace(0, settings["points_time"], settings["frames"]+1)
	freq_xs = settings["cp"] + settings["hub"]*np.sin(2*np.pi*time_xs*settings["fm_freq"])
	fm1_sine_ys = np.sin(1*2*np.pi*settings["fm_freq"]*time_xs + settings["fm1_phase"])
	fm2_sine_ys = np.sin(2*2*np.pi*settings["fm_freq"]*time_xs + settings["fm2_phase"])
	
	# Real Spectrum
	spectrum_xs = np.linspace(-settings["spectrum_width"]/2, settings["spectrum_width"]/2, settings["points_spectrum"])
	spectrum_ys = lineshape(spectrum_xs)
	max_signal = np.nanmax(spectrum_ys)*1.2
	
	# Signal at Detector
	detector_ys = lineshape(freq_xs)
	
	# Resulting Spectra
	result_xs = np.linspace(-settings["spectrum_width"]/2, settings["spectrum_width"]/2, settings["points_spectrum"]+1)
	result_ys_1f = spectrum_dm_modulation(result_xs, lineshape, freq_xs, fm1_sine_ys)*10
	result_ys_2f = spectrum_dm_modulation(result_xs, lineshape, freq_xs, fm2_sine_ys)*10
	
	# Inline FM Sine
	tmp_min = max_signal
	tmp_max = settings["spectrum_height"]
	inline_sine_xs = tmp_min+time_xs*(tmp_max-tmp_min)/settings["points_time"]
	
	# Inline Detector Signal
	tmp_min = settings["inline_sig_start"]
	tmp_max = settings["spectrum_width"]/2
	inline_det_xs = tmp_min+time_xs*(tmp_max-tmp_min)/settings["points_time"]

	# Left Column
	axs[i].plot(spectrum_xs, spectrum_ys, color=color_spectrum, **kwargs_plot)
	axs[i].plot(freq_xs, inline_sine_xs, color=color_fm, linewidth=0.5)
	axs[i].plot(inline_det_xs, detector_ys, color=color_spectrum+"82", linewidth=0.5)
	
	axs[i].set_xlim(-settings["spectrum_width"]/2, settings["spectrum_width"]/2)
	axs[i].set_ylim(0, 3)
	axs[i].set_yticks((0, 1, 2))
	axs[i].set_xticks((-5, 0, +5))
	
	axs[i].text(settings["cp"]+settings["hub"], settings["spectrum_height"]-0.1, f" FM\n Amp: {hub:.2f}", ha="left", va="top", transform=axs[i].transData, size=5, color=color_fm)
	axs[i].text(+5.1, max(detector_ys)+0.05, "Signal at\nDetector", ha="left", va="bottom", transform=axs[i].transData, size=5, color=color_spectrum+"82")
	axs[i].arrow(settings["cp"], max_signal, 0, settings["spectrum_height"]-max_signal-0.1, **arrow_dict)
	axs[i].arrow(4, (max(detector_ys)+min(detector_ys))/2, 5.9, 0, **arrow_dict)
	
	lxs = [settings["cp"]+x*settings["hub"] for x in (0, 0.5, 1)]
	lys = [lineshape(x) for x in lxs]
	axs[i].vlines(lxs, lys, +100, color = "#c2c2c2"+"82", linewidth = 0.35)
	axs[i].hlines(lys, lxs, +100, color = "#c2c2c2"+"82", linewidth = 0.35)
	
	axs[i].text(settings["cp"], settings["spectrum_height"]-0.1, "  $t$  ", ha="right", va="top", transform=axs[i].transData, **fd)
	axs[i].text(9.9, (max(detector_ys)+min(detector_ys))/2-0.1, "$t$ ", ha="right", va="top", transform=axs[i].transData, **fd)
	
	# Middle Column
	axs[10+i].yaxis.tick_right()
	axs[10+i].plot(result_xs, result_ys_1f, color=color_result, **kwargs_plot)
	axs[10+i].set_xticks((-5, 0, +5))
	axs[10+i].set_xlim(-settings["spectrum_width"]/2, settings["spectrum_width"]/2)
	
	axs[10+i].text(1, 0.95, f"Max: {max(result_ys_1f):.2f}  ", ha="right", va="top", transform=axs[10+i].transAxes, size=6)
	axs[10+i].axhline(0, color = "#c2c2c2"+"82", linewidth = 0.5)
	
	# Left Column
	axs[20+i].yaxis.tick_right()
	axs[20+i].plot(result_xs, result_ys_2f, color=color_result_2f, **kwargs_plot)
	axs[20+i].set_xticks((-5, 0, +5))
	axs[20+i].set_xlim(-settings["spectrum_width"]/2, settings["spectrum_width"]/2)
	
	axs[20+i].text(1, 0.95, f"Max: {max(result_ys_2f):.2f}  ", ha="right", va="top", transform=axs[20+i].transAxes, size=6)
	axs[20+i].axhline(0, color = "#c2c2c2"+"82", linewidth = 0.5)
	
	# Set Values for y-Ranges
	if result_ys_1f.min() < res_range_1f[0]:
		res_range_1f[0] = result_ys_1f.min()
	if result_ys_1f.max() > res_range_1f[1]:
		res_range_1f[1] = result_ys_1f.max()
		
	if result_ys_2f.min() < res_range_2f[0]:
		res_range_2f[0] = result_ys_2f.min()
	if result_ys_2f.max() > res_range_2f[1]:
		res_range_2f[1] = result_ys_2f.max()

##
## y-Ranges
##
tmp = res_range_1f[1] - res_range_1f[0]
res_range_1f = [res_range_1f[0]-0.2*tmp, res_range_1f[1]+0.2*tmp]
tmp = res_range_2f[1] - res_range_2f[0]
res_range_2f = [res_range_2f[0]-0.2*tmp, res_range_2f[1]+0.2*tmp]
res_range_1f[1] = res_range_1f[0]*res_range_2f[1]/res_range_2f[0]

for i in range(10, 13):
	axs[i].set_ylim(res_range_1f)
for i in range(20, 23):
	axs[i].set_ylim(res_range_2f)

##
## Column Titles and Axis Labels
##
tmp_fd = {'fontsize': 8}
axs[0].set_title("Real Lineshape", color=color_spectrum, fontdict=tmp_fd)
axs[10].set_title("1f-Lineshape", color=color_result, fontdict=tmp_fd)
axs[20].set_title("2f-Lineshape", color=color_result_2f, fontdict=tmp_fd)
_ = shared_labels(fig, gs[:, :], "Frequency [A.U.]", "Intensity [A.U.]")


##
## Saving/showing
##
fig.subplots_adjust(right = 0.9, top = 0.90, bottom=0.15, left=0.12)
fig.savefig("FMAmplitude.pdf")
# plt.show()
plt.close()