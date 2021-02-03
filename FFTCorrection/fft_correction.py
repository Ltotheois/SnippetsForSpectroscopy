#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Luis Bonah
# Description : FFT Baseline Correction

import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib import gridspec
import scipy.fftpack


from datetime import datetime
from dateutil.parser import parse


def crawl_info_file(filename):
		with open(filename, "r") as file:
			data = file.read()

		blocks = [block.split("\n") for block in data.split("\n\n")]

		db_dict={}
		db_dict["filename"]		=	os.path.basename(filename)
		if db_dict["filename"].find("fat.")!=-1:
			db_dict["status"]	=	0
		elif db_dict["filename"][0]=="B":
			db_dict["status"]	=	+2
		else:
			db_dict["status"]	=	+1
		
		#Block 0
		#General
		db_dict["experimentalist"]	=	blocks[0][0].split(" measured",1)[0]
		db_dict["molecule"]			=	blocks[0][0].split(": ",1)[1]
		[db_dict["datestart"],db_dict["datestop"]]=blocks[0][1].replace("Measurement started on ","").split(" and finished at ")
		db_dict["measurementtype"]	=	blocks[0][2].split(": ",1)[1]
		
		for i in ["datestart", "datestop"]:
			db_dict[i]=str(parse(db_dict[i]))
			
		#Block 1
		#Probe Synthesizer
		db_dict["freqstart"]	=	float(blocks[1][1].split(": ",1)[1].replace("MHz",""))
		db_dict["freqcenter"]	=	float(blocks[1][2].split(": ",1)[1].replace("MHz",""))
		db_dict["freqstop"]		=	float(blocks[1][3].split(": ",1)[1].replace("MHz",""))
		db_dict["span"]			=	float(blocks[1][4].split(": ",1)[1].replace("MHz",""))
		db_dict["steps"]		=	float(blocks[1][5].split(": ",1)[1])
		db_dict["stepsize"]		=	float(blocks[1][6].split(": ",1)[1].replace("kHz",""))
		db_dict["datapoints"]	=	float(blocks[1][7].split(": ",1)[1])
		db_dict["mfactor"]		=	float(blocks[1][8].split(": ",1)[1])
		db_dict["fmdeviation"]	=	float(blocks[1][9].split(": ",1)[1].replace("kHz/V",""))
		db_dict["probepower"]	=	float(blocks[1][10].split(": ",1)[1].replace("dBm",""))

		#Block 2     
		#Pump Synthesizer           
							
		#Block 3                 
		#Lock In                  
		db_dict["delaytime"]	=	float(blocks[3][1].split(": ",1)[1].replace("ms",""))
		db_dict["timeconstant"]	=	float(blocks[3][2].split(": ",1)[1].replace("ms",""))
		db_dict["averagedpoints"]=	float(blocks[3][3].split(": ",1)[1])
		db_dict["averagediter"]	=	float(blocks[3][4].split(": ",1)[1])
		db_dict["oscfreq"]		=	float(blocks[3][5].split(": ",1)[1].replace("Hz",""))
		db_dict["oscamplitude"]	=	float(blocks[3][6].split(": ",1)[1].replace("V",""))
		db_dict["ADCslope"]		=	float(blocks[3][7].split(": ",1)[1].replace("dB/ocatve",""))
		db_dict["ACgain"]		=	float(blocks[3][8].split(": ",1)[1].replace("dB",""))
								
		#Block 3                 
		#Pressure                   
		db_dict["totFMmod"]		=	float(blocks[4][0].split("= ",1)[1].replace("kHz",""))
		if blocks[4][1].split(": ",1)[0]=="Pressure":
			db_dict["pressurestart"]=	"pressure not available"
			db_dict["pressureend"]	=	"pressure not available"
		else:
			db_dict["pressurestart"]=	blocks[4][1].split(": ",1)[1].replace("mbar","")
			db_dict["pressureend"]	=	blocks[4][2].split(": ",1)[1].replace("mbar","")
		for i in ("pressurestart", "pressureend"):
			if db_dict[i].find("pressure not available")!=-1:
				db_dict[i]=-1
			else:
				db_dict[i]=float(db_dict[i])
		
		return(db_dict)


def decrease_standingwave(data, save_data):
	global xs, ys, ys_corr, duration, filename, current_index
	
	
	def onclick(event):
		if event.inaxes == ax1:
			if event.button == 1:
				if isinstance(event.xdata,np.float64):
					cut_off_slider.set_val(event.xdata)
	
	def press(key):
		global xs, ys, ys_corr, duration, filename, current_index
		
		if key=="left":
			cut_off_slider.set_val(cut_off_slider.val-0.01)
		elif key=="shift+left":
			cut_off_slider.set_val(cut_off_slider.val-0.1)
		elif key=="right":
			cut_off_slider.set_val(cut_off_slider.val+0.01)
		elif key=="shift+right":
			cut_off_slider.set_val(cut_off_slider.val+0.1)
		elif key=="up":
			cut_off_slider.set_val(cut_off_slider.val+0.15)
		elif key=="shift+up":
			cut_off_slider.set_val(cut_off_slider.val+0.2)
		elif key=="down":
			cut_off_slider.set_val(cut_off_slider.val-0.15)
		elif key=="shift+down":
			cut_off_slider.set_val(cut_off_slider.val-0.2)
		
		elif key in [" ", "space", "enter"]:
			save_data(xs, ys_corr, filename)
			current_index += 1
			if current_index >= len(data):
				current_index = len(data)-1
			update_plot()
		elif key in ["ctrl+left"]:
			current_index -= 1
			if current_index < 0:
				current_index = 0
			update_plot()
		elif key in ["ctrl+q"]:
			fig.canvas.mpl_disconnect(cid_1)
			fig.canvas.mpl_disconnect(cid_2)
			plt.close()
		elif key in ["escape"]:
			current_index += 1
			if current_index >= len(data):
				current_index = len(data)-1
			update_plot()
	
	def update_plot(rescale = True):
		global xs, ys, ys_corr, duration, filename, current_index
		
		xs, ys, duration, filename = data[current_index]
		cutoff_freq = cut_off_slider.val
		
		fft_ys		= scipy.fftpack.rfft(ys)
		fft_xs		= scipy.fftpack.rfftfreq(len(ys), duration/len(xs))
		
		fft_cut = [x for x in fft_ys]
		fft_bas = [x for x in fft_ys]

		fft_cut = [fft_ys[i] if fft_xs[i] > cutoff_freq else 0 for i in range(len(fft_ys))]
		fft_bas = [fft_ys[i] if fft_xs[i] < cutoff_freq else 0 for i in range(len(fft_ys))]
				
		ys_corr = scipy.fftpack.irfft(fft_cut)
		ys_base = scipy.fftpack.irfft(fft_bas)
		
		ax1.lines[0].set_data(fft_xs, fft_ys)
		ax2.lines[0].set_data(xs, ys)
		ax2.lines[1].set_data(xs, ys_base)
		ax3.lines[0].set_data(xs, ys_corr)
		
		if rescale == True:
			for ax in (ax1, ax2, ax3):
				ax.relim()
				ax.autoscale_view(True, True, True)
		
			tmp = (np.nanmin(fft_xs), np.nanmax(fft_xs))
			tmp_p = (tmp[1]-tmp[0])*0.05
			
			cut_off_slider.valmin = tmp[0]-tmp_p
			cut_off_slider.valmax = tmp[1]+tmp_p
			cut_off_slider.ax.set_xlim(cut_off_slider.valmin, cut_off_slider.valmax)
		
		else:
			ax3.relim()
			ax3.autoscale_view(True, True, True)
		
		line.set_xdata(cutoff_freq)
		fig.suptitle(f"{current_index+1}/{len(data)}: {os.path.basename(filename)}", transform=fig.transFigure, ha="center")

		fig.canvas.draw_idle()
		
	
	current_index = 0
	cutoff_freq = 0.0
	
	fig= plt.figure()
	gs = gridspec.GridSpec(7, 4, height_ratios = [0.5,1,.5,1,1,0.5,0.5], hspace = 0, wspace=0) 

	ax0 = fig.add_subplot(gs[0, :])
	cut_off_slider = Slider(ax0, "Cut-Off", 0, 1, valinit=cutoff_freq)
	cut_off_slider.on_changed(lambda a: update_plot(rescale=False))

	ax1 = fig.add_subplot(gs[1, :])
	ax1.plot([], [], color="green", label="FFT Coefficients")
	ax1.legend(loc = "upper right")
	line = ax1.axvline(x=cutoff_freq, color="red", ls="--")
	
	tmp_ax = fig.add_subplot(gs[2, :])
	tmp_ax.axis("off")

	ax2 = fig.add_subplot(gs[3, :])
	ax2.plot([], [], color="#6ebeff", label="Original Spectrum")	
	ax2.plot([], [], color="#FF0266", label="Baseline", linewidth=3, alpha=0.3)
	ax2.get_xaxis().set_visible(False)
	ax2.legend(loc = "upper right")
	
	ax3 = fig.add_subplot(gs[4, :], sharex=ax2)
	ax3.plot([], [], color="#0336FF", label="Corrected Spectrum")
	ax3.legend(loc = "upper right")
	
	tmp_ax = fig.add_subplot(gs[5, :])
	tmp_ax.axis("off")
	
	buttons = [("Quit", "ctrl+q"), ("Previous", "ctrl+left"), ("Next", "escape"), ("Save", "enter")]
	refs = {}
	
	for i, (text, key) in enumerate(buttons):		
		tmp_ax = fig.add_subplot(gs[6, i])
		tmp_button = Button(tmp_ax, text)
		tmp_button.on_clicked(lambda a, key=key: press(key))
		refs[key] = tmp_button
	
	update_plot()

	cid_1 = fig.canvas.mpl_connect('button_press_event', onclick)
	cid_2 = fig.canvas.mpl_connect('key_press_event', lambda event: press(event.key))

	fig.tight_layout()
	plt.show()

	
if __name__ == '__main__':
	
	# Set up what should happen with corrected data, here just save to file
	def save_data(xs, ys, filename):
		fname, extension = os.path.splitext(filename)
		df = pd.DataFrame({"x": xs, "y":ys})
		df.to_csv(fname + "FFT" + extension, header=False, index=False, sep="\t")
	
	# Get files
	try:
		import tkinter as tk
		from tkinter import filedialog
		root = tk.Tk()
		root.withdraw()
		filenames = filedialog.askopenfilename(multiple=True) 
		root.destroy()
	except Exception as E:
		filenames = input("Enter filenames: ").split(",")
	
	filenames = list(set(filenames))
	if len(filenames) != 0:
	
		# Fill data array
		data = []
		
		for filename in filenames:
			# Get x- and y-data
			df = pd.read_csv(filename, sep="\t", skip_blank_lines=True, dtype=np.float64, names=(["x", "y"]))
			xs = df["x"].to_numpy()
			ys = df["y"].to_numpy()
			
			# Get duration if possible, otherwise set to 30
			fname, extension = os.path.splitext(filename)
			try:
				db_dict=crawl_info_file(fname+".info")
				date_start=parse(db_dict["datestart"])
				date_stop=parse(db_dict["datestop"])
				duration=(date_stop-date_start).total_seconds()/2
			except Exception as E:
				duration = 30
				
			data.append((xs, ys, duration, filename))
		
		# Start main function
		decrease_standingwave(data, save_data)