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
	"figsize" :						(30/2.54, 20/2.54),
	"dpi" :							300,
}
kwargs_plot = {
	"linewidth" :					1,
}

kwargs_font = {
	'size' :						7,
}
matplotlib.rc('font', **kwargs_font)

settings = {
	"splitting" 		:	4,
	"width_gauss"		:	0.5,
	"width_lorentz"		:	0.5,
	"amp"				:	np.float64(1/(2*0.20870928052036772)),
	"fm_freq"			:	np.float64(0.06),
	"am_freq"			:	np.float64(0.01),
	"fm_phase"			:	np.pi/2,
	"am_phase"			:	0,
	"hub"				:	1.5,
	"cp"				:	1.0,
	"duration"			:	100,
	"frames"			:	1000,
	"inline_sig_start"	:	4,
	"inline_fm_start"	:	1.2,
	"pixels_spectrum"	:	1000,
	"spectrum_width"	:	20,
	"spectrum_height"	:	2,
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
	ax.set_xlabel(xlabel, labelpad = 5)
	ax.set_ylabel(ylabel, labelpad = 10)
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
		ind = len(xs)
	if ind <= 0:
		segs = (((0,0),(0, 0)) for i in range(0))
	else:
		segs = (((xs[i], ys[i]),(xs[i+1], ys[i+1])) for i in range(ind))
	coll = matplotlib.collections.LineCollection(segs, color=colors, **kwargs_plot)
	ref = ax.add_collection(coll)
	return(ref)

def spectrum_2f(offsets, lineshape, freq_xs, fm2_sine_ys):
	out = []
	for offset in offsets:
		tmp = lineshape(freq_xs+offset-settings["cp"])*fm2_sine_ys
		res = np.sum(tmp)/len(tmp)
		out.append(res)
	return(np.array(out))

##
## Basic Setup
##

labels = [1, "h1", 2, 3, 4, 5, 6, 7, "h2", 8]
ratios = [1 if str(labels[i]).isnumeric() else 0.6 for i in range(len(labels))]
ratios[3] = 0.5
ratios[6] = 0.5

gs = matplotlib.gridspec.GridSpec(len(labels), 1, height_ratios = ratios, hspace=0)

fig = plt.figure(**kwargs_figure)
axs = {}
for x, label in enumerate(labels):
	axs[label] = fig.add_subplot(gs[x,:])

for x in [1,2,3,4,5,6,8]:
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
mixed2_ys = lockin_ys*am_sine_ys

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

_ = shared_labels(fig, gs[2:8, :], "Time [A.U.]", "Intensity [A.U.]")


# --- First Row ---

# Inline FM Sine
axs[1].plot(freq_xs, inline_sine_xs, color=color_help, linewidth=0.3)

# Inline Signal
plot_inlinedet = make_collection(axs[1], [], [], [], kwargs_plot)

# Spectrum
plot_spectrum, = axs[1].plot(spectrum_xs, spectrum_ys_off, color=color_off, **kwargs_plot)

# Highlight Points and Lines
point_fm,   = axs[1].plot(0, 0, marker='o', markersize=3, color=color_help+"82", zorder=10)
point_spec, = axs[1].plot(0, 0, marker='o', markersize=3, color=color_off,       zorder=10)
point_det,  = axs[1].plot(0, 0, marker='o', markersize=3, color=color_off+"82",  zorder=10)
line_mm1,   = axs[1].plot((0, 0), (0, 0), color=color_help+"82", linewidth=0.5)
line_mm2,   = axs[1].plot((0, 0), (0, 0), color=color_help+"82", linewidth=0.5)
line_cp,    = axs[1].plot((0, 0), (0, 0), color=color_help, zorder=10, ls=":", lw=0.5)

# Styling
axs[1].set_xlim(-settings["spectrum_width"]/2, settings["spectrum_width"]/2)
axs[1].set_ylim(0, settings["spectrum_height"])
axs[1].set_xlabel("Frequency [A.U.]")
axs[1].set_ylabel("Int. [A.U.]", labelpad = 10)
axs[1].text(1.0, 0.5, "  1: Spectrum", ha="left", va="center", transform=axs[1].transAxes)
axs[1].text(+settings["hub"], 1.95, " FM ", ha="left", va="top", transform=axs[1].transData, size=8, color="#c2c2c2")
label_det = axs[1].text(+4, max_signal*1.2, "Signal at Detector", ha="left", va="center", transform=axs[1].transData, size=8, color=color_off+"82")


# --- Hidden Row ---

axs["h1"].axis("off")


# --- Second Row ---

# Detector Signal
plot_detector = make_collection(axs[2], [], [], [], kwargs_plot)

# Styling 
axs[2].set_xlim(0, settings["duration"])
axs[2].set_ylim(0, max_signal)
axs[2].text(1.0, 0.5, "  2: Signal at Detector $D(t)$", ha="left", va="center", transform=axs[2].transAxes)
axs[2].axhline(y=0, color=color_0, linewidth = 0.25)


# --- Third Row ---

# 2f Signal
axs[3].plot(time_xs, fm2_sine_ys, color=color_sin, **kwargs_plot)

# Styling
axs[3].get_yaxis().set_visible(False)
axs[3].set_xlim(0, settings["duration"])
axs[3].set_ylim((-1.8, 1.8))
axs[3].axhline(y=0, color=color_0, linewidth = 0.25)
axs[3].text(1.0, 0.5, r"  3: $\sin(2\cdot 2\pi\cdot f_{FM}\cdot t + \phi)$", ha="left", va="center", transform=axs[3].transAxes)


# --- Fourth Row ---

# Mixed Signal
plot_mixed = make_collection(axs[4], [], [], [], kwargs_plot)

# Filling
zero_arr = np.array([0])
fill1 = axs[4].fill_between(zero_arr, zero_arr, zero_arr, zero_arr>0, color=color_pos)
fill2 = axs[4].fill_between(zero_arr, zero_arr, zero_arr, zero_arr<0, color=color_neg)

# Styling
axs[4].set_xlim(0, settings["duration"])
axs[4].set_ylim(np.nanmin(mixed_ys)*1.2 if np.nanmin(mixed_ys) < 0 else -1, np.nanmax(mixed_ys)*1.2 if np.nanmax(mixed_ys) > 0 else 1)
axs[4].axhline(y=0, color='#878787', linewidth = 0.25)
axs[4].text(1.0, 0.5, r"  4: $\sin(2 \cdot 2\pi\cdot f_{FM}\cdot t + \phi)\cdot D(t)$", ha="left", va="center", transform=axs[4].transAxes)


# --- Fifth Row ---

# LockIn Signal
plot_lockin = make_collection(axs[5], [], [], [], kwargs_plot)

# Integration Points
points_int = axs[5].scatter([], [], color=[], s=12)

# Styling
axs[5].set_xlim(0, settings["duration"])
axs[5].set_ylim(np.nanmin(lockin_ys)*1.4 if np.nanmin(lockin_ys) < 0 else -1, np.nanmax(lockin_ys)*1.4 if np.nanmax(lockin_ys) > 0 else 1)
axs[5].text(1.0, 0.5, r"  5: Output of Lock-In", ha="left", va="center", transform=axs[5].transAxes)

# Arrows between plots
kwargs_int_arrows_1 = {
	"axesA"		: axs[4],
	"axesB"		: axs[5],
	"coordsA"	: "data",
	"coordsB"	: "data",
	"arrowstyle": "-",
	"color"		: "#878787",
}

kwargs_int_arrows_2 = {
	"axesA"		: axs[7],
	"axesB"		: axs[8],
	"coordsA"	: "data",
	"coordsB"	: "data",
	"arrowstyle": "-",
	"color"		: "#878787",
}


# --- Sixth Row ---

# AM Signal
axs[6].plot(time_xs, am_sine_ys, color=color_sin, **kwargs_plot)

# Styling
axs[6].get_yaxis().set_visible(False)
axs[6].set_xlim(0, settings["duration"])
axs[6].set_ylim((-1.8, 1.8))
axs[6].axhline(y=0, color=color_0, linewidth = 0.25)
axs[6].text(1.0, 0.5, r"  6: $\sin(2\pi\cdot f_{AM}\cdot t + \varphi)$", ha="left", va="center", transform=axs[6].transAxes)


# --- Seventh Row ---

# Mixed 2 Lockin Signal
plot_mixed2 = make_collection(axs[7], [], [], [], kwargs_plot)

# Filling
zero_arr = np.array([0])
fill3 = axs[7].fill_between(zero_arr, zero_arr, zero_arr, zero_arr>0, color=color_pos)
fill4 = axs[7].fill_between(zero_arr, zero_arr, zero_arr, zero_arr<0, color=color_neg)

# Styling
axs[7].set_xlim(0, settings["duration"])
axs[7].set_ylim(np.nanmin(mixed2_ys)*1.2 if np.nanmin(mixed2_ys) < 0 else 0, np.nanmax(mixed2_ys)*1.2 if np.nanmax(mixed2_ys) > 0 else 1)
axs[7].axhline(y=0, color=color_0, linewidth = 0.25)
axs[7].text(1.0, 0.5, r"  7: $\sin(2\pi\cdot f_{AM}\cdot t + \varphi)\cdot S(t)$", ha="left", va="center", transform=axs[7].transAxes)


# --- Hidden Row ---

axs["h2"].axis("off")


# --- Eigth Row ---

# DMDR Signal
tmp_xs = np.linspace(-settings["spectrum_width"]/2, settings["spectrum_width"]/2, 101)
tmp_ys = -spectrum_2f(tmp_xs, lineshape_on, freq_xs, fm2_sine_ys)
axs[8].plot(tmp_xs, tmp_ys, color=color_dmdr+"40", **kwargs_plot)
plot_dmdr = axs[8].scatter([], [], color=color_dmdr, s=12)

# Styling
axs[8].set_xlim(-settings["spectrum_width"]/2, settings["spectrum_width"]/2)
axs[8].set_ylim(np.nanmin(tmp_ys)*1.2, np.nanmax(tmp_ys)*1.2)
axs[8].set_xlabel("Frequency [A.U.]")
axs[8].set_ylabel("Int. [A.U.]", labelpad = 10)
axs[8].text(1.0, 0.5, "  8: DMDR Spectrum", ha="left", va="center", transform=axs[8].transAxes)


##
## Animation
##

last_integration = 0
last_integration_i = 0
integration_points = []
integration_colors = []
flashing_artists = []
def animate(i):
	ext_upd = []
	if i<settings["frames"]:
		ind = i+1
		if am_sine_ys[i] > 0:
			spectrum_ys = spectrum_ys_off
			color_ooo = color_off
		else:
			spectrum_ys = spectrum_ys_on
			color_ooo = color_on
		
		x = freq_xs[i]
		y = detector_ys[i]
		
		global last_integration, last_integration_i, integration_colors, integration_points
		global flashing_artists
		global plot_inlinedet, plot_spectrum, label_det, point_fm, point_det, point_spec, line_mm1, line_mm2, plot_detector, plot_mixed, fill1, fill2, plot_lockin, plot_mixed2, fill3, fill4, points_int, line_cp
		
		# First Row
		plot_inlinedet.remove()
		plot_inlinedet = make_collection(axs[1], inline_det_xs, detector_ys, [x+"82" for x in colors], kwargs_plot, ind=ind)
		
		plot_spectrum.set_data(spectrum_xs, spectrum_ys)
		plot_spectrum.set_color(color_ooo)
		label_det.set_color(color_ooo+"82")
		
		point_fm.set_data(x, inline_sine_xs[i])
		point_det.set_data(inline_det_xs[i], y)
		point_det.set_color(color_ooo+"82")
		point_spec.set_data(x, y)
		point_spec.set_color(color_ooo+"82")
		
		line_cp.set_data((settings["cp"], settings["cp"]), (0, 100))
		line_mm1.set_data([x, x], [y, 100])
		line_mm2.set_data([x, 100], [y, y])
		
		# Second Row
		plot_detector.remove()
		plot_detector = make_collection(axs[2], time_xs, detector_ys, colors, kwargs_plot, ind=ind)
		
		# Fourth Row
		plot_mixed.remove()
		plot_mixed = make_collection(axs[4], time_xs, mixed_ys, colors, kwargs_plot, ind=ind)
		
		fill1.remove()
		fill2.remove()
		
		fill1 = axs[4].fill_between(time_xs[:ind], 0*time_xs[:ind], mixed_ys[:ind], mixed_ys[:ind]>0, color=color_pos)
		fill2 = axs[4].fill_between(time_xs[:ind], 0*time_xs[:ind], mixed_ys[:ind], mixed_ys[:ind]<0, color=color_neg)
		
		i_i = time_xs[i]// (1/settings["fm_freq"])
		if i_i != last_integration or i == settings["frames"]-1:
			
			# Fifth Row
			plot_lockin.remove()
			plot_lockin = make_collection(axs[5], time_xs, lockin_ys, colors, kwargs_plot, ind=ind)
		
			# Seventh Row
			plot_mixed2.remove()
			plot_mixed2 = make_collection(axs[7], time_xs, mixed2_ys, colors, kwargs_plot, ind=ind)
			
			fill3.remove()
			fill4.remove()
			
			fill3 = axs[7].fill_between(time_xs[:ind], 0*time_xs[:ind], mixed2_ys[:ind], mixed2_ys[:ind]>0, color=color_pos)
			fill4 = axs[7].fill_between(time_xs[:ind], 0*time_xs[:ind], mixed2_ys[:ind], mixed2_ys[:ind]<0, color=color_neg)
			
			# Arrows
			flashing_artists = []
			
			last_x = time_xs[last_integration_i+1]
			this_x = time_xs[i]
			mid_x = (last_x + this_x)/2
			mid_y = lockin_ys[int((last_integration_i+i)/2)]
			
			flashing_artists.append(axs[4].axvspan(last_x, this_x, facecolor=color_help+"30", alpha=0.5))
			
			patch1 = matplotlib.patches.ConnectionPatch(xyA=(last_x, np.nanmin(mixed_ys)*1.2), xyB=(mid_x, mid_y), **kwargs_int_arrows_1)
			patch2 = matplotlib.patches.ConnectionPatch(xyA=(this_x, np.nanmin(mixed_ys)*1.2), xyB=(mid_x, mid_y), **kwargs_int_arrows_1)
			axs[5].add_artist(patch1)
			axs[5].add_artist(patch2)
			flashing_artists.append(patch1)
			flashing_artists.append(patch2)
			
			integration_points.append((mid_x, mid_y))
			integration_colors.append(colors[int((last_integration_i+i)/2)])
			points_int.set_offsets(integration_points)
			points_int.set_color(integration_colors)
			
			ext_upd.extend(flashing_artists)
			
			# clean up 
			last_integration = i_i
			last_integration_i = i
		elif (i - last_integration_i) < 20 and last_integration_i != 0:
			ext_upd.extend(flashing_artists)
		else:
			for artist in flashing_artists:
				try:
					artist.remove()
				except:
					pass
		if i == settings["frames"]-1:
			ext_upd.append(axs[7].axvspan(0, settings["duration"], facecolor=color_help+"30", alpha=0.5))
			
			patch3 = matplotlib.patches.ConnectionPatch(xyA=(0, np.nanmin(mixed2_ys)*1.2), xyB=resulting_point, **kwargs_int_arrows_2)
			patch4 = matplotlib.patches.ConnectionPatch(xyA=(settings["duration"], np.nanmin(mixed2_ys)*1.2), xyB=resulting_point, **kwargs_int_arrows_2)
			plot_dmdr.set_offsets(resulting_point)
			ext_upd.append(plot_dmdr)
			axs[8].add_artist(patch3)
			axs[8].add_artist(patch4)
			ext_upd.append(patch3)
			ext_upd.append(patch4)
		
		# Artists to be updated
	
	update_patches = [plot_inlinedet, plot_spectrum, label_det, point_fm, point_det, point_spec, line_mm1, line_mm2, plot_detector, plot_mixed, fill1, fill2, plot_lockin, plot_mixed2, fill3, fill4, points_int, line_cp] + ext_upd
	
	return update_patches


##
## Main Work
##

save = False
show = False
if len(sys.argv) > 1:
	save = True
	interval = 30
else:
	show = True
	interval = 0


fig.subplots_adjust(right = 0.8, top = 0.96)
anim = animation.FuncAnimation(fig, animate, frames=settings["frames"]+100, interval=interval, blit=True, repeat=False, save_count=settings["frames"]+100)
if save:
	anim.save("animation_dmdr.mp4", writer="ffmpeg")
if show:
	plt.show()
plt.close()