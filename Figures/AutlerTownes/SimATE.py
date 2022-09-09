#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Luis Bonah

import numpy as np
import scipy.special
import matplotlib
import matplotlib.pyplot as plt

c = {
	"pump_frequency": 50050,
	"xrange": (76000, 76100),
	"xpoints": 1000,
	"sigma": 1,
	"gamma": 0.2,
	"fm_frequency": 27000,
	"fm_amplitude": 0.2,
	"fm_timeconstant": 0.050,
	"fm_points": 1000,
	"efield": 2,
	"color_off": "#6ebeff",
	"color_on": "#FF0266",
	"color_dmdr": "#000000",
}

fig_kwargs = {
	"figsize": (2.5, 1.5),
	"dpi": 300,
}

def lineshape(xs, x0, y0, sigma, gamma):
	return (scipy.special.voigt_profile(xs-x0, sigma, gamma)*y0)

def calculate_signals(peaks):
	xs = np.linspace(*c["xrange"], c["xpoints"])
	
	# Real Spectrum
	ys = np.sum([lineshape(xs, x0, y0, c["sigma"], c["gamma"]) for x0, y0, *_ in peaks], axis=0)
	
	# Dr Spectrum
	dr_ys = np.sum([lineshape(xs, x0, y0, c["sigma"], c["gamma"]) for x0, y0, *_ in dr_peaks(peaks)], axis=0)
	
	sine_xs = np.linspace(0, c["fm_timeconstant"], c["fm_points"])
	sine_ys = np.sin(2*np.pi*sine_xs*c["fm_frequency"])*c["fm_amplitude"]
	
	# Demodulated Spectrum
	demod_ys = []
	for x in xs:
		tmp = np.interp(sine_ys + x, xs, ys)*np.cos(2*np.pi*sine_xs*2*c["fm_frequency"])
		demod_ys.append(np.mean(tmp))
	
	# Demodulated DR Spectrm
	demod_dr_ys= []
	for x in xs:
		tmp = np.interp(sine_ys + x, xs, dr_ys)*np.cos(2*np.pi*sine_xs*2*c["fm_frequency"])
		demod_dr_ys.append(np.mean(tmp))
		
	return(xs, ys, dr_ys, sine_xs, sine_ys, np.array(demod_ys), np.array(demod_dr_ys)) 

def dr_peaks(peaks):
	peaks_dr = []
	
	for probe_x, probe_y, pump_x, pump_y, regressive in peaks:
		detuning = c["pump_frequency"] - pump_x
		rabi = c["efield"]*pump_y
		offset_1 = 1/2*(+np.sqrt(detuning**2 + rabi**2)-detuning)
		offset_2 = 1/2*(-np.sqrt(detuning**2 + rabi**2)-detuning)
		
		if regressive:
			offset_1, offset_2 = -offset_1, -offset_2
		
		intfactor = 1 if detuning == 0 else np.tan(0.5 * np.arctan(-rabi/detuning))**2
		int_1, int_2 = 1/(intfactor+1), intfactor/(intfactor+1)
		
		if detuning < 0:
			int_1, int_2 = int_2, int_1
		
		peaks_dr.append((probe_x+offset_1, probe_y*int_1, pump_x, pump_y))
		peaks_dr.append((probe_x+offset_2, probe_y*int_2, pump_x, pump_y))
		
	return(peaks_dr)


def plot(peaks, interactive=False):
	if interactive:
		gs = matplotlib.gridspec.GridSpec(6, 2, wspace=0.1, height_ratios=(1,1,1,0.25,0.25,0.25))
	else:
		gs = matplotlib.gridspec.GridSpec(3, 2, wspace=0.1)

	fig = plt.figure()
	
	tmp = {
		1: gs[0, 0],
		2: gs[0, 1],
		3: gs[1, 0],
		4: gs[1, 1],
		5: gs[2, :],
	}
	axs = {key: fig.add_subplot(pos) for key, pos in tmp.items()}
	for ax in axs.values():
		ax.set_xlim(*c["xrange"])
		ax.set_xticks([])
	
	axs[2].yaxis.tick_right()
	axs[4].yaxis.tick_right()
	
	axs[1].get_shared_x_axes().join(*(ax for ax in axs.values()))
	axs[1].get_shared_y_axes().join(axs[1], axs[2])
	axs[3].get_shared_y_axes().join(axs[3], axs[4])
	
	xs, ys, dr_ys, sine_xs, sine_ys, demod_ys, demod_dr_ys = calculate_signals(peaks)
	
	# Real Spectrum
	axs[1].set_title("Classic")
	axs[1].plot(xs, ys, color=c["color_off"])
	
	# Dr Spectrum
	axs[2].set_title("DR")
	axs[2].plot(xs, dr_ys, color=c["color_on"])
	
	# Demodulated Spectrum
	axs[3].plot(xs, demod_ys, color=c["color_off"])
	
	# Demodulated DR Spectrm
	axs[4].plot(xs, demod_dr_ys, color=c["color_on"])
	
	# DMDR Spectrum
	axs[5].plot(xs, demod_ys-demod_dr_ys, color=c["color_dmdr"])
	
	if interactive:
		tmpax1 = fig.add_subplot(gs[3, :])
		tb1 = matplotlib.widgets.Slider(tmpax1, "Distance", 0, 10, valinit=dist)
		tb1.on_changed(lambda x: update("dist", x))
		
		tmpax2 = fig.add_subplot(gs[4, :])
		tb2 = matplotlib.widgets.Slider(tmpax2, "Fm Amp", 0, 10, valinit=c["fm_amplitude"])
		tb2.on_changed(lambda x: update("fmamp", x))
		
		tmpax3 = fig.add_subplot(gs[5, :])
		tb3 = matplotlib.widgets.Slider(tmpax3, "E Field", 0, 10, valinit=c["efield"])
		tb3.on_changed(lambda x: update("efield", x))
		
		def autoscaley(ax, *yss):
			ys = np.concatenate(yss)
			min, max = np.min(ys), np.max(ys)
			dist = max-min
			
			ax.set_ylim(min-0.1*dist, max+0.1*dist)
		
		def update(key, value):
			global c, peaks
			
			if key == "dist":
				peaks = (
					#(probe x, probe y, pump x, pump dipole, regressive)
					(76050-value, 1, 50050-value, 1, False), 
					(76050+value, 1, 50050+value, 1, False),
				)
			elif key == "fmamp":
				c["fm_amplitude"] = value
			elif key == "efield":
				c["efield"] = value
			
			xs, ys, dr_ys, sine_xs, sine_ys, demod_ys, demod_dr_ys = calculate_signals(peaks)

			# Real Spectra
			axs[1].get_lines()[0].set_data(xs, ys)		
			axs[2].get_lines()[0].set_data(xs, dr_ys)
			autoscaley(axs[1], dr_ys, ys)
			autoscaley(axs[2], dr_ys, ys)
			
			# Demodulated Spectra
			axs[3].get_lines()[0].set_data(xs, demod_ys)
			axs[4].get_lines()[0].set_data(xs, demod_dr_ys)
			autoscaley(axs[3], demod_ys, demod_dr_ys)
			autoscaley(axs[4], demod_ys, demod_dr_ys)
			
			# DMDR Spectrum
			axs[5].get_lines()[0].set_data(xs, demod_ys-demod_dr_ys)
			autoscaley(axs[5], demod_ys-demod_dr_ys)
			
		plt.show()
	
	return (fig)

if __name__ == "__main__":
	interactive = True
	
	if interactive:
		dist = 1.1
		peaks = (
			#(probe x, probe y, pump x, pump dipole, regressive)
			(76050-dist, 1, 50050-dist, 1, False), 
			(76050+dist, 1, 50050+dist, 1, False),
		)
		fig = plot(peaks, interactive=True)
	else:
		for fm_amplitude in np.linspace(0, 1, 11):
			c["fm_amplitude"] = fm_amplitude
			for dist in np.linspace(0, 3, 31):
				peaks = (
					#(probe x, probe y, pump x, pump dipole, regressive)
					(76050-dist, 1, 50050-dist, 1, False), 
					(76050+dist, 1, 50050+dist, 1, False),
				)
				fig = plot(peaks)
				fig.savefig(f"fmamp_{fm_amplitude:.1f}_dist_{dist:.1f}.png")
				plt.close(fig)