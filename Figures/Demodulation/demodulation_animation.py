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
from matplotlib.gridspec import GridSpec
from scipy import special


##
## Settings
##
pump_on = False

kwargs_figure = {
	"figsize" :						(30/2.54, 14/2.54),
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
	"splitting" 		:	4,
	"width_gauss"		:	0.5,
	"width_lorentz"		:	0.5,
	"amp"				:	np.float64(1/(2*0.20870928052036772)),
	"fm_freq"			:	np.float64(25)/2/np.pi,
	"am_freq"			:	np.float64(10e-3),
	"hub"				:	1.5,
	"cp"				:	0.0,
	"DR"				:	True,
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
	"size"	: 8,
}

color_sin = "#0336FF"
color_up = "#6ebeff"
color_at = "#FF0266"


##
## Helper Function
##

def shared_labels(fig, gsp, xlabel, ylabel):
	ax = fig.add_subplot(gsp, frameon=False)
	ax.spines['top'].set_color('none')
	ax.spines['bottom'].set_color('none')
	ax.spines['left'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.tick_params(axis="both", labelcolor='#fff0', top=False, bottom=False, left=False, right=False, zorder=-100)
	ax.set_xlabel(xlabel, labelpad = 10)
	ax.set_ylabel(ylabel, labelpad = 20)
	return(ax)

def lineshape(xs):
	ys = special.voigt_profile(xs, settings["width_gauss"], settings["width_lorentz"])
	return(ys*settings["amp"])
	
def lineshape_at(xs):
	ys = 0.5*(special.voigt_profile(xs-settings["splitting"]/2, settings["width_gauss"], settings["width_lorentz"])+special.voigt_profile(xs+settings["splitting"]/2, settings["width_gauss"], settings["width_lorentz"]))
	return(ys*settings["amp"])

##
## Basic Setup
##

nop = 6
nof = 1000
gs = GridSpec(nop, 1, height_ratios = [1, 1, 1, 0.5, 1, 1], hspace=0)

fig = plt.figure(**kwargs_figure)
axs = {}
for x in range(nop):
	axs[x] = fig.add_subplot(gs[x,:])

for x in range(1, nop):
	axs[x].get_xaxis().set_visible(False)

# axs[2].get_shared_x_axes().join(axs[2], axs[3], axs[4], axs[5], axs[6], axs[7], axs[8], axs[9])


##
## Prepare Data
##
sig_xs = np.linspace(-10,10,1000)
time_sin_xs = np.linspace(1.2, 2, nof+1)
time_sin_ys = settings["cp"]+settings["hub"]*np.sin(2*np.pi*settings["fm_freq"]*(time_sin_xs-time_sin_xs[0]))
tmp_xs = np.linspace(4, 10, nof+1)
mul_ys_pre = np.concatenate((lineshape(time_sin_ys)*settings["hub"]*np.sin(np.pi/2+4*np.pi*settings["fm_freq"]*(time_sin_xs-time_sin_xs[0])), lineshape_at(time_sin_ys)*settings["hub"]*np.sin(np.pi/2+4*np.pi*settings["fm_freq"]*(time_sin_xs-time_sin_xs[0]))))

last_integration = 0
last_integration_i = -np.inf
integral = []
integral_ind = set()

i = 0
det_xs = []
det_ys = []
time_xs = []
colors = []

axs[i].plot(time_sin_ys, time_sin_xs, color="#c2c2c2", **kwargs_plot)
segs = (((0,0),(0, 0)) for i in range(0))
coll = matplotlib.collections.LineCollection(segs, color=colors, **kwargs_plot)
det_signal = axs[i].add_collection(coll)
sig_signal, = axs[i].plot(sig_xs, lineshape(sig_xs), color=color_up, **kwargs_plot)
axs[i].set_xlim(-10, 10)
axs[i].set_ylim(0, 2)
axs[i].set_xlabel("Frequency [A.U.]")
axs[i].set_ylabel("Int. [A.U.]", labelpad = 10)
axs[i].text(1.0, 0.5, "  1: Spectrum", ha="left", va="center", transform=axs[i].transAxes)


axs[i].text(+settings["hub"], 1.95, " FM ", ha="left", va="top", transform=axs[i].transData, size=8, color="#c2c2c2")
det_label = axs[i].text(+4, 1.2, "Signal at Detector", ha="left", va="center", transform=axs[i].transData, size=8, color=color_up+"82")

time_point, = axs[i].plot(0, 0,  marker='o', markersize=3, color="#c2c2c2"+"82", zorder=10)
data_point, = axs[i].plot(0, 0,  marker='o', markersize=3, color=color_up, zorder=10)
spec_point, = axs[i].plot(0, 0, marker='o', markersize=3, color=color_up+"82", zorder=10)
mm_line,    = axs[i].plot((0, 0), (0, 0), color="#c2c2c2"+"82", linewidth=0.5)
mm_line2,   = axs[i].plot((0, 0), (0, 0), color="#c2c2c2"+"82", linewidth=0.5)

i += 1
axs[i].axis("off")

i += 1
xs = np.linspace(0,nof, nof+1)
segs = (((0,0),(0, 0)) for i in range(0))
coll = matplotlib.collections.LineCollection(segs, color=colors, **kwargs_plot)
time_plot = axs[i].add_collection(coll)
axs[i].set_xlim(0, nof)
axs[i].set_ylim(0, 1.2)
axs[i].text(1.0, 0.5, "  2: Signal at Detector $D(t)$", ha="left", va="center", transform=axs[i].transAxes)
axs[i].axhline(y=0, color='#878787', linewidth = 0.25)
axs[i+1].set_xlim(0, nof)
axs[i+1].set_ylim((-1.8, 1.8))
axs[i+1].plot(xs, settings["hub"]*np.sin(np.pi/2+2*2*np.pi*settings["fm_freq"]*(time_sin_xs-time_sin_xs[0])), color=color_sin, **kwargs_plot)
axs[i+1].axhline(y=0, color='#878787', linewidth = 0.25)
axs[i+1].text(1.0, 0.5, r"  3: $\sin(2\cdot 2\pi\cdot f_{FM}\cdot t + \phi)$", ha="left", va="center", transform=axs[i+1].transAxes)

i += 2

mul_ys = []
segs = (((0,0),(0, 0)) for i in range(0))
coll = matplotlib.collections.LineCollection(segs, color=colors, **kwargs_plot)
mul_signal = axs[i].add_collection(coll)
axs[i].set_xlim(0, nof)
axs[i].set_ylim(np.amin(mul_ys_pre)*1.2, np.amax(mul_ys_pre)*1.2)
zero_arr = np.array([0])
fill1 = axs[i].fill_between(zero_arr, zero_arr, zero_arr, zero_arr>0, color="#CCCCCC")
fill2 = axs[i].fill_between(zero_arr, zero_arr, zero_arr, zero_arr<0, color="#89919A")
axs[i].axhline(y=0, color='#878787', linewidth = 0.25)
axs[i].text(1.0, 0.5, r"  4: $\sin(2 \cdot 2\pi\cdot f_{FM}\cdot t + \phi)\cdot D(t)$", ha="left", va="center", transform=axs[i].transAxes)

## Arrows between plots
conn_kwargs = {
	"axesA"		: axs[4],
	"axesB"		: axs[5],
	"coordsA"	: "data",
	"coordsB"	: "data",
	"arrowstyle": "-",
	"color"		: "#878787",
}
patch1 = None
patch2 = None

i += 1
output_xs = []
output_ys = []
axs[i].set_xlim(0, nof)
axs[i].set_ylim(-0.5, 0.5)
axs[i].text(1.0, 0.5, r"  5: Output of Lock-In", ha="left", va="center", transform=axs[i].transAxes)
out_signal, = axs[i].plot([], [], color="black", marker='o', markersize=3, linewidth=0)

_ = shared_labels(fig, gs[2:6, :], "Time [A.U.]", "Intensity [A.U.]")


def animate(i):
	global integral
	global last_integration
	global last_integration_i
	
	global lineshape
	global lineshape_at
	global pump_on
	global time_plot
	global det_signal
	global mul_signal
	
	global patch1
	global patch2
	
	if pump_on:
		colors.append(color_at)
	else:
		colors.append(color_up)
	
	x = time_sin_ys[i]
	y = lineshape(x)
	y_mul = y*settings["hub"]*np.sin(np.pi/2+4*np.pi*settings["fm_freq"]*(time_sin_xs[i]-time_sin_xs[0]))
	
	det_xs.append(4+i*np.float(6)/nof)
	det_ys.append(y)
	time_xs.append(i)
	mul_ys.append(y_mul)
	
	mm_line.set_data([x, x], [y, 100])
	mm_line2.set_data([x, 100], [y, y])
	
	time_point.set_data(x, time_sin_xs[i])
	data_point.set_data(x, y)
	spec_point.set_data(tmp_xs[i], y)
	
	det_signal.remove()
	segs = (((det_xs[i],det_ys[i]),(det_xs[i+1],det_ys[1+i])) for i in range(len(det_xs)-1))
	coll = matplotlib.collections.LineCollection(segs, colors=[color+"82" for color in colors])
	det_signal = axs[0].add_collection(coll) 

	time_plot.remove()
	segs = (((time_xs[i],det_ys[i]),(time_xs[i+1],det_ys[1+i])) for i in range(len(time_xs)-1))
	coll = matplotlib.collections.LineCollection(segs, colors=colors)
	time_plot = axs[2].add_collection(coll) 

	
	axs[4].collections.clear()
	segs = (((time_xs[i],mul_ys[i]),(time_xs[i+1],mul_ys[1+i])) for i in range(len(time_xs)-1))
	coll = matplotlib.collections.LineCollection(segs, colors=colors)
	mul_signal = axs[4].add_collection(coll) 

	ys_numpy = np.array(mul_ys)
	xs_numpy = np.array(time_xs)

	fill1 = axs[4].fill_between(xs_numpy, 0*xs_numpy, ys_numpy, ys_numpy>0, color="#CCCCCC")
	fill2 = axs[4].fill_between(xs_numpy, 0*xs_numpy, ys_numpy, ys_numpy<0, color="#89919A")
	tbu = [sig_signal, mm_line, mm_line2, time_point, data_point, spec_point, det_signal, time_plot, mul_signal, fill1, fill2, out_signal]
	
	if i not in integral_ind:
		integral_ind.add(i)
		integral.append(y_mul)
	
	tmp = int((time_sin_xs[i]-time_sin_xs[0])*settings["fm_freq"]/1)
	if tmp > last_integration:
		
		if settings["DR"] and tmp%1==0:
			pump_on = not pump_on
			if pump_on:
				tmp_color = color_at
			else:
				tmp_color = color_up
			sig_signal.set_color(tmp_color)
			det_label.set_color(tmp_color+"82")
			data_point.set_color(tmp_color)
			spec_point.set_color(tmp_color+"82")
			
			lineshape, lineshape_at = lineshape_at, lineshape
			sig_signal.set_data(sig_xs, lineshape(sig_xs))
		res = sum(integral)/len(integral)
		integral = []
		last_integration = tmp
		last_integration_i = i
		output_xs.append(i*(tmp-0.5)/tmp)
		output_ys.append(res)
		# output_ys.extend([res]*(len(time_xs)-len(output_ys)))
		out_signal.set_data(output_xs, output_ys)
		tmp_range = min(output_ys), max(output_ys)
		diff = tmp_range[1]-tmp_range[0]
		
		patch1 = matplotlib.patches.ConnectionPatch(xyA=(i*(tmp-1)/tmp, 0), xyB=(i*(tmp-0.5)/tmp, res), **conn_kwargs)
		patch2 = matplotlib.patches.ConnectionPatch(xyA=(i*(tmp-0)/tmp, 0), xyB=(i*(tmp-0.5)/tmp, res), **conn_kwargs)
		axs[5].add_artist(patch1)
		axs[5].add_artist(patch2)
		tbu.append(patch1)
		tbu.append(patch2)
	elif (i - last_integration_i) < 20:
		tbu.append(patch1)
		tbu.append(patch2)
		
		
	return tbu
	

fig.subplots_adjust(right = 0.8, top = 0.96)
anim = animation.FuncAnimation(fig, animate, frames=nof, interval=0, blit=True, repeat=False, save_count=nof)
anim.save("animation.mp4", writer="ffmpeg")
plt.show()
plt.close()