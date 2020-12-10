#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Luis Bonah
# Description : Loomis-Wood Plot software

import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib import gridspec
import scipy.fftpack

from datetime import datetime
from dateutil.parser import parse

global size
size=(6,6)
return_code = 0

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


def decrease_standingwave(xs, ys, duration, label = "FFT Correction"):
	global ys_corr
	global size
	global return_code
	
	def update_plot(value=None):
		global ys_corr
		cutoff_freq = cut_off_slider.val
		fft_cut = [x for x in fft_ys]
		fft_bas = [x for x in fft_ys]
		for i in range(0,len(fft_cut)):
			if fft_xs[i]<cutoff_freq:
				fft_cut[i] = 0
			else:
				fft_bas[i] = 0
				
		ys_corr = scipy.fftpack.irfft(fft_cut)
		ys_base = scipy.fftpack.irfft(fft_bas)
		ax3.lines[0].set_data(xs, ys_base)
		ax4.lines[0].set_data(xs, ys_corr)
		
		ax3.relim()
		ax4.relim()
		ax3.autoscale_view(True,True,True)
		ax4.autoscale_view(True,True,True)
		line.set_xdata(cutoff_freq)
		fig.canvas.draw_idle()
	
	def onclick(event):
		if event.inaxes==fig.axes[0]:
			if event.button == 1:
				if isinstance(event.xdata,np.float64):
					cut_off_slider.set_val(event.xdata)
			
	def press(event):
		global return_code
		global size
		if event.key=="left":
			cut_off_slider.set_val(cut_off_slider.val-0.01)
		elif event.key=="shift+left":
			cut_off_slider.set_val(cut_off_slider.val-0.1)
		elif event.key=="right":
			cut_off_slider.set_val(cut_off_slider.val+0.01)
		elif event.key=="shift+right":
			cut_off_slider.set_val(cut_off_slider.val+0.1)
		elif event.key=="up":
			cut_off_slider.set_val(cut_off_slider.val+0.15)
		elif event.key=="shift+up":
			cut_off_slider.set_val(cut_off_slider.val+0.2)
		elif event.key=="down":
			cut_off_slider.set_val(cut_off_slider.val-0.15)
		elif event.key=="shift+down":
			cut_off_slider.set_val(cut_off_slider.val-0.2)
		elif event.key in ["ctrl+right", "ctrl+left", "ctrl+q", "escape", "enter", "space", " "]:
			fig.canvas.mpl_disconnect(cid_1)
			fig.canvas.mpl_disconnect(cid_2)
			size=fig.get_size_inches()
			plt.close()
			if event.key in [" ", "space", "enter"]:
				return_code = 1
			elif event.key in ["ctrl+left"]:
				return_code = "back"
			elif event.key in ["ctrl+q"]:
				return_code = -1
			else:
				return_code = 0
	
	initial_cut_off = 0.2
	fft_ys		= scipy.fftpack.rfft(ys)
	fft_xs		= scipy.fftpack.rfftfreq(len(ys), duration/len(xs))
	
	fig= plt.figure(figsize=size)
	gs = gridspec.GridSpec(7, 1, height_ratios = [1,.5,1,1,1,0.5,1], hspace = 0) 


	ax1 = fig.add_subplot(gs[0])
	ax1.plot(fft_xs, fft_ys, color="green", label="FFT Coefficients")
	ax1.legend(loc = "upper right")
	line = ax1.axvline(x=initial_cut_off, color="red", ls="--")
	
	tmp_ax = fig.add_subplot(gs[1])
	tmp_ax.axis("off")

	ax2 = fig.add_subplot(gs[2])
	ax2.plot(xs, ys, color="blue", label="Original Spectrum")
	ax2.get_xaxis().set_visible(False)
	ax2.legend(loc = "upper right")
	
	ax3 = fig.add_subplot(gs[3], sharex=ax2)
	ax3.plot([], [], color="green", label="Baseline")
	ax3.get_xaxis().set_visible(False)
	ax3.legend(loc = "upper right")
	
	ax4 = fig.add_subplot(gs[4], sharex=ax2)
	ax4.plot([], [], color="red", label="Corrected Spectrum")
	ax4.legend(loc = "upper right")
	
	tmp_ax = fig.add_subplot(gs[5])
	tmp_ax.axis("off")

	ax6= fig.add_subplot(gs[6])
	cut_off_slider = Slider(ax6, "Cut-Off Freq", *ax1.get_xlim(), valinit=initial_cut_off)
	cut_off_slider.on_changed(update_plot)
	
	update_plot()

	cid_1 = fig.canvas.mpl_connect('button_press_event', onclick)
	cid_2 = fig.canvas.mpl_connect('key_press_event', press)
	fig.suptitle(os.path.basename(label), transform=fig.transFigure, ha="center")
	fig.tight_layout()
	plt.show()
	
	return(xs, ys_corr, return_code)
	
if __name__ == '__main__':
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
	
	i = 0
	while i < len(filenames):
		filename = filenames[i]
		df = pd.read_csv(filename, sep="\t", skip_blank_lines=True, dtype=np.float64, names=(["x", "y"]))
		xs = df["x"].to_numpy()
		ys = df["y"].to_numpy()
		
		fname, extension = os.path.splitext(filename)
		try:
			db_dict=crawl_info_file(fname+".info")
			date_start=parse(db_dict["datestart"])
			date_stop=parse(db_dict["datestop"])
			duration=(date_stop-date_start).total_seconds()/2
		except Exception as E:
			duration = 30
		xs, ys, rc = decrease_standingwave(xs, ys, duration, filename)
		if return_code == 1:
			df = pd.DataFrame({"x": xs, "y":ys})
			df.to_csv(fname + "FFT" + extension, header=False, index=False, sep="\t")
		elif return_code == -1:
			break
		elif return_code == "back":
			i -= 2
		i += 1
