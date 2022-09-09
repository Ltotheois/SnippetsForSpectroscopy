#!/usr/bin/env python 
#-*- coding: utf-8 -*-

import os
import traceback
import numpy as np
import plotly.offline as py
import plotly.graph_objs as go

def get_pump_frequency_from_filename(file):
	return(float(file.split("_")[3]))

def get_filenames(folder=None, fileending=".dat"):
	if folder == None:
		import tkinter as tk
		from tkinter import filedialog
		root = tk.Tk()
		root.withdraw()
		folder = filedialog.askdirectory() 
		root.destroy()
	filenames = [file for file in os.listdir(folder) if file.endswith(fileending)]
	return(filenames,folder)

def get_data_from_files(filenames, folder):
	data_dict={}
	for file in filenames:
		try:
			pump_frequency = get_pump_frequency_from_filename(file)
		except ValueError:
			print(f"Could not find the pump-frequency for {file}.")
			continue
		tmp_data = np.genfromtxt(os.path.join(folder, file))
		tmp_data = sorted(tmp_data,key= lambda line:line[0])
		data_dict[pump_frequency] = [line[1] for line in tmp_data]
	y_values = sorted(data_dict.keys())
	x_values = sorted([line[0] for line in tmp_data])
	z_matrix = []
	for i in sorted(data_dict.keys(), reverse = False):
		z_matrix.append(data_dict[i])
	return(x_values, y_values, z_matrix)

def plot_surface(x_values, y_values, z_matrix, filename=None):
	filename = "2Dspectroscopy.html" if filename==None else filename

	layout = go.Layout(
		title="",
		scene = dict(
						xaxis = dict(
							title="",
							ticktext= x_values[0::50],
							tickvals= x_values[0::50]),
						yaxis = dict(
							title="",
							ticktext= y_values[0::50],
							tickvals= y_values[0::50]),
						zaxis = dict(
							title="")),
				autosize=True,
				width=1200,
				height=800,
				margin=dict(
						l=65,
						r=50,
						b=65,
						t=90
				)
	)
	
	updatemenus=list([
		dict(
			buttons=list([
				dict(
					args=['type', 'surface'],
					label='3D Surface',
					method='restyle'
				),
				dict(
					args=['type', 'heatmap'],
					label='Heatmap',
					method='restyle'
				)
			]),
			direction = 'left',
			pad = {'r': 10, 't': 10},
			showactive = True,
			type = 'buttons',
			x = 0.1,
			xanchor = 'left',
			y = 1.1,
			yanchor = 'top'
		),
	])
	layout['updatemenus'] = updatemenus

	# Only Z-matrix and no x and y values will improve the animation while interacting
	#fig = go.Figure(data=[go.Surface(z=z_matrix)], layout=layout)
	
	# Including x and y will give you the proper freqs for each point, but make animations less enjoyable
	fig = go.Figure(data=[go.Surface(z=z_matrix, x=x_values, y=y_values)], layout=layout)
	
	py.plot(fig, filename=filename)


try:
	if __name__=="__main__":
		plot_surface(*get_data_from_files(*get_filenames()))
except Exception as E:
	print("An error occurred!")
	print(E)
	print()
	print(traceback.format_exc())
	input()
