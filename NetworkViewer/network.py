#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Luis Bonah
# Description: Show Network and According Measurements

import sys, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms

## Parameters
##

config = {
	"acolor": "#0336FF",
	"bcolor": "#FF0266",
	"awcolor": "#6ebeff",
	"canvas_xrange": (+0.0,  9.0),
	"canvas_yrange": (-0.5, 12.5),
	"separator": "\t",
	"offset": 0.1,
}
kwargs_plot = {
	"linewidth" :	1.5,
	"linestyle" :	"solid",
	"color" :		"#000000",
}
kwargs_figure = {
	"figsize" :		(15/2.54, 10/2.54),
	"dpi" :			300,
}
kwargs_font = {
	# 'family' :	'normal',
	# 'weight' :	'bold',
	'size' :		10,
}

arrow_dict = {
	"head_width" :		0.1,
	"head_length" :		0.15,
	"length_includes_head" : True,
	"width" :			0.01,
	"picker":			10,
}

matplotlib.rc('font', **kwargs_font)

dict_annotation = {
	"horizontalalignment": 			"right",
	"verticalalignment": 			"center",
	"fontdict" :					{"size" : 4},
}

states = [
	(( 8,0, 8),	1.3, 0.0),
	(( 9,0, 9),	1.3, 1.0),
	((10,0,10),	1.3, 2.5),
	((11,0,11),	1.3, 4.6),
	((12,0,12),	1.3, 7.4),
	((13,0,13),	1.3,10.4),
	
	(( 8,1, 8),	3.4, 0.4),
	(( 9,1, 9),	3.4, 1.4),
	((10,1,10),	3.4, 2.9),
	((11,1,11),	3.4, 5.0),
	((12,1,12),	3.4, 7.8),
	((13,1,13),	3.4,10.8),
	
	(( 8,1, 7),	6.2, 1.1),
	(( 9,1, 8),	6.2, 2.1),
	((10,1, 9),	6.2, 3.6),
	((11,1,10),	6.2, 5.7),
	((12,1,11),	6.2, 8.5),
	((13,1,12),	6.2,11.5),
	
	(( 8,2, 7),	8.3, 1.5),
	(( 9,2, 8),	8.3, 2.5),
	((10,2, 9),	8.3, 4.0),
	((11,2,10),	8.3, 6.1),
	((12,2,11),	8.3, 8.9),
	((13,2,12),	8.3,11.9),
]

transitions = [
	(( 9,0, 9), ( 8,0, 8), "qR0/DR_88057.98_Pump_79430.5+79310.5_18-Jan-21-6.26.49 PM.dat"),
	((10,0,10), ( 9,0, 9), "qR0/DR_88057.98_Pump_79430.5+79310.5_18-Jan-21-6.26.49 PM.dat"),
	((11,0,11), (10,0,10), "qR0/DR_96635.44_Pump_88057.98+87937.98_18-Jan-21-6.28.36 PM.dat"),
	((12,0,12), (11,0,11), "qR0/DR_105166.675_Pump_96635.44+96515.44_18-Jan-21-6.30.27 PM.dat"),
	((13,0,13), (12,0,12), "qR0/DR_113657.775_Pump_105166.675+105046.675_18-Jan-21-6.31.21 PM.dat"),
	
	(( 9,1, 9), ( 8,1, 8), "qR1/DR_78000_Pump_92319.2832+92199.2832_09-Apr-21-1.31.51 PM_onlyForward.dat"),
	((10,1,10), ( 9,1, 9), "qR1/DR_87000_Pump_78103.9+77983.9_09-Apr-21-3.14.41 PM_onlyForward.dat"),
	((11,1,11), (10,1,10), "qR1/DR_95271_Pump_86697.7253+86577.7253_09-Apr-21-3.27.02 PM_onlyForward.dat"),
	((12,1,12), (11,1,11), "qR1/DR_104000_Pump_95285.3785+95165.3785_09-Apr-21-3.59.47 PM_onlyForward.dat"),
	((13,1,13), (12,1,12), "qR1/DR_112800_Pump_103864.1459+103744.1459_09-Apr-21-4.27.10 PM_onlyForward.dat"),

	((10,0,10), ( 9,1, 9), "pR1/DR_L1_75169.2_Pump_86697.7253+86577.7253_09-Apr-21-4.41.38 PM_onlyForward.dat"),
	((11,0,11), (10,1,10), "pR1/DR_L2_85106.9147_Pump_95285.3785+95165.3785_09-Apr-21-4.42.23 PM_onlyForward.dat"),
	((12,0,12), (11,1,11), "pR1/DR_L3_94988.2112_Pump_103864.1459+103744.1459_09-Apr-21-4.43.07 PM_onlyForward.dat"),
	((13,0,13), (12,1,12), "pR1/DR_L4_104781.8403_Pump_112432.4613+112312.4613_09-Apr-21-4.43.52 PM_onlyForward.dat"),
	
	(( 9,1, 9), ( 8,0, 8), "rR0/DR_92500_Pump_79430.5+79310.5_09-Apr-21-10.31.31 AM_onlyForward.dat"),
	((10,1,10), ( 9,0, 9), "rR0/DR_L1_99586.5053_Pump_88057.98+87937.98_09-Apr-21-4.35.30 PM_onlyForward.dat"),
	((11,1,11), (10,0,10), "rR0/DR_L2_106813.9038_Pump_96635.44+96515.44_09-Apr-21-4.37.15 PM_onlyForward.dat"),
	((12,1,12), (11,0,11), "rR0/DR_L3_114042.6097_Pump_105166.675+105046.675_09-Apr-21-4.39.00 PM_onlyForward.dat"),
	((13,1,13), (12,0,12), "rR0/DR_L4_121308.396_Pump_113657.775+113537.775_09-Apr-21-4.40.46 PM_onlyForward.dat"),
]

## Plotting figure
##
def main():
	fig, ax = plt.subplots(1, 1, **kwargs_figure)
	ax.get_xaxis().set_visible(False)
	ax.get_yaxis().set_visible(False)
	ax.axis("off")
	ax.axis('equal')
	
	ax.set_xlim(*config["canvas_xrange"])
	ax.set_ylim(*config["canvas_yrange"])
	
	for state in states:
		xs = np.linspace(-0.5, +0.5, 101) + state[1]
		ys = xs*0+state[2]
		ax.plot(xs, ys, **kwargs_plot)
		label = f"${state[0][0]}_{{{state[0][1]}, {state[0][2]}}}$"
		ax.text(**dict_annotation, x=state[1]-0.6, y=state[2]+0.05, s=label)
	
	def get_coords(state):
		for tmp in states:
			if state == tmp[0]:
				return(tmp)

	def get_color(transition):
		delta_J = transition[0][0] - transition[1][0]
		delta_Ka = transition[0][1] - transition[1][1]
		delta_Kc = transition[0][2] - transition[1][2]
		
		if delta_J == 1 and delta_Ka == 0 and delta_Kc != 0:
			return(config["acolor"])
		elif delta_J == 0 and delta_Ka == 0 and delta_Kc != 0:
			return(config["awcolor"])
		else:
			return(config["bcolor"])
	
	def get_offset(start, end):
		if start > end:
			return -config["offset"]
		elif end > start:
			return +config["offset"]
		else:
			return 0
	
	arrows = {}
	for transition in transitions:
		_, x1, y1 = get_coords(transition[0])
		_, x2, y2 = get_coords(transition[1])
		color = get_color(transition)
		
		xoffset = get_offset(x2, x1)
		yoffset = get_offset(y2, y1)
		
		arrow = ax.arrow(x2+xoffset, y2+yoffset, x1-x2-2*xoffset, y1-y2-2*yoffset, color=color, **arrow_dict)
		arrows[arrow] = transition
	
	def onclick(event):
		transition = arrows[event.artist]
		fig, ax = plt.subplots()
		data = np.genfromtxt(transition[2], delimiter=config["separator"])
		xs = data[:, 0]
		ys = data[:, 1]
		ax.plot(xs, ys)
		fig.show()

	cid = fig.canvas.mpl_connect('pick_event', onclick)
	
	fig.tight_layout()
	plt.show()
	
if __name__ == "__main__":
	main()