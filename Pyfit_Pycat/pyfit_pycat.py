#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Luis Bonah
# Description : SPFIT/SPCAT wrapping Library


import os
import argparse
import subprocess
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

## Miscellaneous
preset_colors = ['#000000e0', '#dbb243e0', '#2e42d3e0', '#e54fe3e0', '#f23434e0']

## Cat Dataframe Format
_cat_df_columns = ["x", "error", "y", "degfreed", "elower", "usd", "tag", "qnfmt", 'qnu1', 'qnu2', 'qnu3', 'qnu4', 'qnu5', 'qnu6', 'qnl1', 'qnl2', 'qnl3', 'qnl4', 'qnl5', 'qnl6', 'comment']
_cat_df_dtypes = [np.float64, np.float64, np.float64, np.int16, np.float64, np.int16, np.int32, np.int16]+[np.int16]*12+[str]
_cat_widths = [13, 8, 8, 2, 10, 3, 7, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 100]

## Lin Dataframe Format
_lin_df_columns = ['qnu1', 'qnu2', 'qnu3', 'qnu4', 'qnu5', 'qnu6', 'qnl1', 'qnl2', 'qnl3', 'qnl4', 'qnl5', 'qnl6', "x", "error", "weight", "comment"]
_lin_df_dtypes = [np.int16]*12 + [np.float64, np.float64, np.float64]+[str]
_lin_widths = range(0,37,3)

## Egy Dataframe Format
_egy_df_columns = ['iblk', 'indx', 'egy', 'err', 'pmix', 'we', ':', 'qn1', 'qn2', 'qn3', 'qn4', 'qn5', 'qn6']
_egy_df_dtypes = [np.int16, np.int16, np.float64, np.float64, np.float64, np.int16, str]+[np.int16]*6
_egy_df_widths = [0,6,11,29,47,58,63,64,67,70,73,76,79,82]

## Helper
def column_to_float(val):
	val = val.strip()
	if val == "" or val == ":":
		return(-1)
	elif val[0].isalpha():
		val = str(ord(val[0].upper())-55)+val[1:]
	return(np.float64(val))

def tmp_dir_name(dir):
	i = 0
	tmp_dir = os.path.join(dir, f"pyfit_pycat_temp_{i}")
	while os.path.isdir(tmp_dir):
		i += 1
		tmp_dir = os.path.join(dir, f"pyfit_pycat_temp_{i}")
	return(tmp_dir)

def rmdir(dir):
	for file in os.listdir(dir):
		os.remove(os.path.join(dir, file))
	os.rmdir(dir)

## Format Functions
def lin_to_df(fname):
	dtypes = _lin_df_dtypes
	widths = _lin_widths
	column_names = _lin_df_columns
	dtypes_dict = {column_names[i]:dtypes[i] for i in range(len(column_names))}
	
	data = []
	with open(fname, "r") as file:
		for line in file:
			if line.strip() == "" or line.startswith("#"):
				continue
			tmp_line_content = [column_to_float(line[i:j]) for i,j in zip(widths[:-1], widths[1:])]+[column_to_float(x) for x in line[36:].split()]
			tmp_line_content.extend([-1]*(15-len(tmp_line_content)))
			data.append(tmp_line_content)
	data = pd.DataFrame(data, dtype = np.float64)
	data["filename"] = fname
	data.columns = column_names
	data = data.astype(dtypes_dict)
	
	## Set correct columns for data
	qn_labels = column_names[0:12]
	QNs = len(qn_labels)
	for i in range(len(qn_labels)):
		tmp_unique = data[qn_labels[i]].unique()
		if len(tmp_unique) == 1 and tmp_unique[0] == -1:
			QNs = i
			break
	QNs = int(QNs/2)
	columns_qn = [f"qnu{i+1}" for i in range(QNs)]+[f"qnl{i+1}" for i in range(QNs)]+[f"qnu{i+1}" for i in range(QNs, 6)]+[f"qnl{i+1}" for i in range(QNs, 6)]
	new_column_names = columns_qn + column_names[12:]
	data.columns = new_column_names

	data.sort_values("x", inplace = True)
	return(data)

def cat_to_df(fname):
	dtypes = _cat_df_dtypes
	widths = _cat_widths
	column_names = _cat_df_columns
	dtypes_dict = {column: dtype for column, dtype in zip(column_names, dtypes)}
	converters = {key: column_to_float for key in column_names[8:20]}
	for key in column_names[8:20]:
		del dtypes_dict[key]
	
	data = pd.read_fwf(fname, widths=widths, names=column_names, converters=converters, skip_blank_lines=True, dtype=dtypes_dict)
	data["filename"] = fname
	data["y"] = 10 ** data["y"]
	
		
	data.sort_values("x", inplace = True)
	return(data)

def df_to_cat(df):
	output = []

	for index, row in df.iterrows():
		freq = row["x"]
		if freq > 99999999.9999:
			freq = 99999999.9999
		elif freq < 0.0001:
			freq = 0.0001
		error = row["error"]
		if error > 999.9999:
			error = 999.9999
		elif error < -99.9999:
			error = -99.9999
		intens = np.log10(row["y"]) if row["y"] > 0 else 0
		if intens < -99.999:
			intens = -99.999
		elif intens > 999.999:
			intens = 999.999

		QNs_string = ""
		for QNlabel in ['qnu1', 'qnu2', 'qnu3', 'qnu4', 'qnu5', 'qnu6', 'qnl1', 'qnl2', 'qnl3', 'qnl4', 'qnl5', 'qnl6']:
			if row[QNlabel] in ["", -1]:
				QNs_string += "  "
			else:
				QNs_string += f"{row[QNlabel]:2.0f}"

		output.append(f"{freq:13.4f}{error:8.4f}{intens:8.4f}{row['degfreed']:2.0f}{row['elower']:10.4f}{row['usd']:3.0f}{row['tag']:7.0f}{row['qnfmt']:4.0f}{QNs_string}{r'  // '+row['comment'] if row['comment'] not in ['', -1, np.NaN] else ''}")
	
	output = "\n".join(output)
	return(output)

def df_to_lin(df):
	output = []
	for index, row in df.iterrows():
		QNs_string = ""
		pad_string = ""
		for QNlabel in ['qnu1', 'qnu2', 'qnu3', 'qnu4', 'qnu5', 'qnu6', 'qnl1', 'qnl2', 'qnl3', 'qnl4', 'qnl5', 'qnl6']:
			if row[QNlabel] in ["", -1]:
				pad_string += "   "
			else:
				QNs_string += f"{row[QNlabel]:3.0f}"
		QNs_string = QNs_string + pad_string
		output.append(f"{QNs_string} {row['x']:13.4f} {row['error']:8.4f} {row['weight']:13.4f}{r'  // '+row['comment'] if row['comment'] not in ['', -1, np.NaN] else ''}")
	
	output = "\n".join(output)
	return(output)

def egy_to_df(fname):
	dtypes = _egy_df_dtypes
	catindices = _egy_df_widths
	column_names = _egy_df_columns
	dtypes_dict = {column_names[i]:dtypes[i] for i in range(len(column_names))}

	data = []
	with open(fname, "r") as file:
		for line in file:
			if line.strip() == "" or line.startswith("#"):
				continue
			data.append([column_to_float(line[i:j]) for i,j in zip(catindices[:-1], catindices[1:])])
	data = pd.DataFrame(data, dtype = np.float64)

	data.columns = column_names
	data = data.astype(dtypes_dict)
	return(data)

def parvar_to_dict(fname):
	result = {}
	with open(fname, "r") as file:
		result["TITLE"] = file.readline().replace("\n", "")
		
		keys = ['NPAR', 'NLINE', 'NITR', 'NXPAR', 'THRESH ', 'ERRTST', 'FRAC', 'CAL']
		result.update({key: value for key, value in zip(keys, file.readline().split())})
		
		keys = ['CHR', 'SPIND', 'NVIB', 'KNMIN', 'KNMAX', 'IXX', 'IAX', 'WTPL', 'WTMN', 'VSYM', 'EWT', 'DIAG', 'XOPT']
		result.update({key: value for key, value in zip(keys, file.readline().split())})
		
		for key, value in result.items():
			if key not in ["TITLE", "CHR"]:
				value = np.float64(value)
				if value%1 == 0:
					result[key] = int(value)
				else:
					result[key] = value
		
		result['STATES'] = []
		if result['VSYM'] < 0:
			for x in range(abs(result['NVIB'])-1):
				line = file.readline()[1:]
				keys = ['SPIND', 'NVIB', 'KNMIN', 'KNMAX', 'IXX', 'IAX', 'WTPL', 'WTMN', 'VSYM', 'EWT', 'DIAG', 'XOPT'] #Only their in case list is changed to dict
				stateline = [int(value) for key, value in zip(keys, line.split())]
				result['STATES'].append(stateline)
				if stateline[8] > 0:
					break
		
		result['PARAMS'] = []
		for line in file:
			keys = ["IDPAR", "PAR", "ERPAR", "LABEL"] #Only their in case list is changed to dict
			funcs = [int, np.float64, np.float64, lambda x: x.replace("/", "")]
			paramline = [func(value) for key, value, func in zip(keys, line.split(), funcs)]
			
			result['PARAMS'].append(paramline)
			
	return(result)

def dict_to_parvar(dct):
	output = []
	output.append(dct["TITLE"])
	
	formats = ['{:4.0f}', ' {:7.0f}', ' {:5.0f}', ' {:4.0f}', '   {: .4e}', '   {: .4e}', '   {: .4e}', ' {:13.4f}']
	
	values = [dct[key] for key in ['NPAR', 'NLINE', 'NITR', 'NXPAR', 'THRESH ', 'ERRTST', 'FRAC', 'CAL'] if key in dct]
	line = "".join([fs.format(x) for x, fs in zip(values, formats)])
	output.append(line)
	
	formats = [' {:4.0f}', ' {:3.0f}', ' {:3.0f}', ' {:4.0f}', ' {:4.0f}', ' {:4.0f}', ' {:4.0f}', ' {:4.0f}', ' {: 7.0f}', ' {:4.0f}', ' {:1.0f}', ' {:4.0f}']
	
	values = [dct[key] for key in ['NVIB', 'KNMIN', 'KNMAX', 'IXX', 'IAX', 'WTPL', 'WTMN', 'VSYM', 'EWT', 'DIAG', 'XOPT'] if key in dct ]
	line = f"{dct['SPIND']}"+ "".join([fs.format(x) for x, fs in zip(values, formats)])
	output.append(line)
	
	for state in dct["STATES"]:
		line = "".join([fs.format(x) for x, fs in zip(state, formats)])
		output.append(line)
	
	for param in dct['PARAMS']:
		comment = ""
		if len(param) > 3:
			comment = f"/{param[3]}"
		output.append(f"{param[0]:13} {param[1]: .15e} {param[2]: .8e} {comment}")
	
	output = "\n".join(output)
	return(output)

def int_to_dict(fname):
	result = {}
	with open(fname, "r") as file:
		result["TITLE"] = file.readline().replace("\n", "")
		
		keys = ['FLAGS', 'TAG', 'QROT', 'FBGN', 'FEND', 'STR0', 'STR1', 'FQLIM', 'TEMP', 'MAXV']
		funcs = [int, int, np.float64, int, int, np.float64, np.float64, np.float64, np.float64, int]
		result.update({key: func(value) for key, value, func in zip(keys, file.readline().split(), funcs)})
		
		result['INTS'] = []
		for line in file:
			keys = ['IDIP', 'DIPOLE']
			funcs = [int, np.float64]
			intline = [func(value) for key, value, func in zip(keys, line.split(), funcs)]
			
			result['INTS'].append(intline)
	
	return(result)

def dict_to_int(dct):
	output = []
	output.append(dct["TITLE"])
	
	formats = ['{:4.0f}', ' {:7.0f}', ' {:13.4f}', ' {:4.0f}', ' {:4.0f}', ' {: 6.2f}', ' {: 6.2f}', ' {:13.4f}', ' {:13.4f}', ' {:4.0f}']
	
	values = [dct[key] for key in ['FLAGS', 'TAG', 'QROT', 'FBGN', 'FEND', 'STR0', 'STR1', 'FQLIM', 'TEMP', 'MAXV'] if key in dct]
	line = "".join([fs.format(x) for x, fs in zip(values, formats)])
	output.append(line)
	
	for param in dct['INTS']:
		output.append(f" {param[0]: d}  {param[1]:.2f}")
	
	output = "\n".join(output)
	return(output)

## Helper Functions
def run_spcat(filename, parameterfile="", path="", wd=None):
	command = f"{os.path.join(path, 'spcat')} {filename} {parameterfile}"
	return(run_subprocess(command, wd))

def run_spfit(filename, parameterfile="", path="", wd=None):
	command = f"{os.path.join(path, 'spfit')} {filename} {parameterfile}"
	return(run_subprocess(command, wd))

def run_subprocess(command, wd = None):
	if wd != None:
		cwd = os.getcwd()
		os.chdir(wd)

	output = subprocess.check_output(command)
	output = output.decode("utf-8")
	
	if wd != None:
		os.chdir(wd)
	
	return(output)

def run_spfit_v(par_dict, lin_df, spfit_path):
	cwd = os.getcwd()
	tmp_dir = tmp_dir_name(cwd)
	
	os.mkdir(tmp_dir)
	os.chdir(tmp_dir)
	
	with open("tmp.par", "w+") as par_file, open("tmp.lin", "w+") as lin_file:
		lin_file.write(df_to_lin(lin_df))
		par_file.write(dict_to_parvar(par_dict))
	
	message = run_spfit("tmp", path=spfit_path, wd=tmp_dir)
	
	result = {"message": message}
	for ext in (".bak", ".par", ".var", ".fit", ".bin"):
		if os.path.isfile(f"tmp{ext}"):
			with open(f"tmp{ext}", "r") as file:
				result[ext] = file.read()
	
	os.chdir(cwd)
	rmdir(tmp_dir)
	return(result)

def run_spcat_v(var_dict, int_dict, spcat_path):
	cwd = os.getcwd()
	tmp_dir = tmp_dir_name(cwd)
	
	os.mkdir(tmp_dir)
	os.chdir(tmp_dir)
	
	with open("tmp.int", "w+") as int_file, open("tmp.var", "w+") as var_file:
		int_file.write(dict_to_int(int_dict))
		var_file.write(dict_to_parvar(var_dict))
	
	message = run_spcat("tmp", path=spcat_path, wd=tmp_dir)
	
	result = {"message": message}
	for ext in (".out", ".cat", ".str", ".egy"):
		if os.path.isfile(f"tmp{ext}"):
			with open(f"tmp{ext}", "r") as file:
				result[ext] = file.read()
	
	os.chdir(cwd)
	rmdir(tmp_dir)
	return(result)

def parse_fit_result(msg):
	results = {}
	
	i_0 = msg.rfind("MICROWAVE RMS")
	i_1 = msg.find(",", i_0)
	RMS = msg[i_0:i_1].split("=")[1].split()[0]
	RMS = np.float64(RMS)
	results["RMS"] = RMS
	
	i_0 = msg.rfind("RMS ERROR=")
	i_1 = msg.find("\n", i_0)
	WRMS = msg[i_0:i_1].split()[-1]
	results["WRMS"] = WRMS
	
	return(results)

def plot_bars(ys, xlabels, title="", xlabel="", ylabel=""):
	fig, ax = plt.subplots()
	plt.title(title)
	plt.ylabel(xlabel)
	plt.xlabel(ylabel)
	
	xs = np.arange(len(ys))
	colors = preset_colors[:len(ys)]
	
	plt.bar(xs, ys, align="center", color=colors)
	plt.xticks(rotation=45, ha='right')

	ax.set_xticks(xs)
	ax.set_xticklabels(xlabels)
	
	return(fig)

## Main Actions
def residuals(df_cat, df_lin, query_string=None, scatter_dict={}, head_number=20, save_fname=None):
	df_merge = pd.merge(df_cat, df_lin, how="inner", on=[f"{tag}{i}" for tag in ("qnu", "qnl") for i in range(1, 7)])
	df_merge["diff"] = df_merge["x_x"] - df_merge["x_y"]
	df_merge["absdiff"] = abs(df_merge["diff"])

	if query_string != None:
		df_merge.query(query_string, inplace=True)
	
	xs = df_merge["x_y"].to_numpy()
	ys = df_merge["diff"].to_numpy()
	
	fig, ax = plt.subplots()
	sc = ax.scatter(xs, ys, **scatter_dict)
	ax.set_xlabel("Frequency")
	ax.set_ylabel("Diff")
	
	annot = ax.annotate("", xy=(0,0), xytext=(5, 5), textcoords="offset points", bbox=dict(boxstyle="round", fc="w"), fontsize=8)
	annot.set_visible(False)
	
	def update_annot(ind):
		pos = sc.get_offsets()[ind["ind"][0]]
		annot.xy = pos
		
		text = []
		for i in ind["ind"]:
			row = df_merge.loc[i]
			upper_state = []
			lower_state = []
			for j in range(6):
				qnu = row[f"qnu{j+1}"]
				qnl = row[f"qnl{j+1}"]
				if qnu != -1 and qnl != -1:
					upper_state.append(str(int(qnu)))
					lower_state.append(str(int(qnl)))
			upper_state = ", ".join(upper_state)
			lower_state = ", ".join(lower_state)
			text.append(f"{upper_state} $\leftarrow$ {lower_state}")
		
		text = "\n".join(text)
		annot.set_text(text)


	def hover(event):
		vis = annot.get_visible()
		if event.inaxes == ax:
			cont, ind = sc.contains(event)
			if cont:
				update_annot(ind)
				annot.set_visible(True)
			else:
				annot.set_visible(False)					
		else:
			annot.set_visible(False)
		fig.canvas.draw_idle()
	
	fig.canvas.mpl_connect("motion_notify_event", hover)
	
	
	df_merge.sort_values("absdiff", inplace=True, ascending=False)
	
	print(df_merge[["x_x", "x_y", "diff", 'qnu1', 'qnu2', 'qnu3', 'qnu4', 'qnu5', 'qnu6', 'qnl1', 'qnl2', 'qnl3', 'qnl4', 'qnl5', 'qnl6']].head(head_number))
	
	plt.tight_layout()
	if save_fname == None:
		plt.show()
	elif type(save_fname) == str:
		plt.savefig(save_fname)

	plt.close()

	return(df_merge)

def energy_levels(egy_df, states, kas, Jmax=60, red_fac=(4714.188+4235.085)/2, save_fname=None, annotate=False, kwargs_plot={}, colors_dict={}):
	fig, ax = plt.subplots()
	
	for ka in kas:
		for state in states:
			for qnsum in (0, 1):
				tmp_df = egy_df.query(f"qn2 == {ka} and qn4 == {state} and qn1 + {qnsum} == qn2+qn3 and qn1 < {Jmax}").copy()
				tmp_df["yV"] = (tmp_df["egy"]*29979.2458-tmp_df["qn1"]*(tmp_df["qn1"]+1)*red_fac)/1e6
				xs = tmp_df["qn1"].to_numpy()
				ys = tmp_df["yV"].to_numpy()
				color = colors_dict.get(state, "#000000")
				marker = "." if qnsum == 0 else "+"
				
				ax.plot(xs, ys, color=color, marker = marker, **kwargs_plot)
				if len(xs) > 0 and qnsum==0 and annotate==True:
					ax.text(x=xs[0], y=ys[0], verticalalignment="center", horizontalalignment="right", transform=ax.transData, fontdict={"color":color, "size":8}, s=f"$K_{{a}}={ka}$  ")
	
	ax.set_xlabel(r"$J$")
	ax.set_ylabel(r"$E_{red}$ [$10^{3}$ GHz]")
	
	plt.tight_layout()
	if save_fname == None:
		plt.show()
	elif type(save_fname) == str:
		plt.savefig(path_save)

	plt.close()

def check_crossings(egy_df, states, kas, Jmax=60, save_fname=None):
	output = []
	series_list = []
	for state in states:
		for ka in kas:
			qnsums = (0, 1) if ka != 0 else (0,)
			for qnsum in qnsums:
				tmp_df = egy_df.query(f"qn2 == {ka} and qn4 == {state} and qn1 + {qnsum} == qn2+qn3 and qn1 < {Jmax}").copy()
				xmin = tmp_df["qn1"].to_numpy().min() if len(tmp_df["qn1"].to_numpy()) != 0 else 0
				ys = tmp_df["egy"].to_numpy()
				
				series_list.append(((state, ka, qnsum), xmin, ys))
	
	crossings = []
	for i in range(len(series_list)):
		for j in range(i+1, len(series_list)):
			desc_1, xmin_1, ys_1 = series_list[i]
			desc_2, xmin_2, ys_2 = series_list[j]
			
			xdiff = xmin_1 - xmin_2
			xstart = max(xmin_1, xmin_2)
			
			if xdiff > 0:
				ys_2 = ys_2[xdiff:]
			elif xdiff < 0:
				ys_1 = ys_1[abs(xdiff):]
			
			ydiff = ys_1 - ys_2
			ytmp = [ydiff[k]*ydiff[k+1] for k in range(len(ydiff)-1)]
			
			for k in range(len(ytmp)):
				if ytmp[k] < 0:
					crossings.append((desc_1, desc_2, xstart+k))
	
	crossings = sorted(crossings, key = lambda x: (x[0][0], x[1][0]))
	
	output.append("Format is state1, state2 @ ka1, J-ka1-kc1 & ka2, J-ka2-kc2 @ Ji, Jf")
	for crossing in crossings:
		J = crossing[2]
		output.append(f"{crossing[0][0]:3d}, {crossing[1][0]:3d} @ {crossing[0][1]:3d}, {crossing[0][2]:1d} & {crossing[1][1]:3d}, {crossing[1][2]:1d} @ {J:3d}, {J+1:3d}")
	
	output.append(f"\nFound {len(crossings)} crossings in total.")
	output = "\n".join(output)
	
	if save_fname == None:
		print(output)
	elif type(save_fname) == str:
		with open(save_fname, "w+") as file:
			file.write(output)

def mixing_coefficient(egy_df, query_string, save_fname=None):
	gs = matplotlib.gridspec.GridSpec(1, 3, width_ratios = [1,0.2, 0.1], hspace=0, wspace=0)
	fig = plt.figure()

	ax = fig.add_subplot(gs[0,0])
	eax = fig.add_subplot(gs[0,1])
	eax.axis("off")
	cbaxs = fig.add_subplot(gs[0,2])
	
	ax.set_xlabel("$J$")
	ax.set_ylabel("$K_{a}$")
	
	tmp_df = egy_df.query(query_string).copy()
	xs = tmp_df["qn1"].to_numpy()
	ys = tmp_df["qn2"].to_numpy()
	zs = tmp_df["pmix"].abs().to_numpy()
	
	df = pd.DataFrame({"x": xs, "y": ys, "z": zs})
	zmatrix = df.pivot_table(values="z", index="y", columns="x")
	zmatrix = zmatrix.to_numpy()

	if len(xs) == 0 or len(ys) == 0:
		print("No data found.")
		return()

	xs = [x-0.5 for x in sorted(list(set(xs)))]
	xs.append(max(xs)+1)
	ys = [y-0.5 for y in sorted(list(set(ys)))]
	ys.append(max(ys)+1)

	clim = (0.5,1)
	ax.pcolormesh(xs, ys, zmatrix)
	ax.set_xlim(min(xs), max(xs))
	ax.set_ylim(min(xs), max(ys))
	
	norm = matplotlib.colors.Normalize(vmin=0.5,vmax=1)
	sm = plt.cm.ScalarMappable(cmap="plasma_r", norm=norm)
	sm.set_array([])
	cb = fig.colorbar(sm, cax=cbaxs, orientation="vertical")
	cb.set_label('Mixing Coefficient', labelpad=10)
	
	plt.tight_layout()
	if save_fname == None:
		plt.show()
	elif type(save_fname) == str:
		plt.savefig(save_fname)
	
	plt.close()

def add_parameter(par_dict, lin_df, param_candidates, spfit_path=None, save_fname=None):
	runs = []
	orig_par_dict = par_dict.copy()
	
	param_candidates.insert(0, ["Initial", 0, 0])
	
	for i, param in enumerate(param_candidates):
		par_dict =  orig_par_dict.copy()
		if param[0] != "Initial":
			par_dict["PARAMS"].append(param)
		
		results = run_spfit_v(par_dict, lin_df, spfit_path)
		RMS = parse_fit_result(results["message"])["RMS"]
		runs.append([i, RMS, par_dict["PARAMS"].copy()])
	
	runs = sorted(runs, key=lambda x: x[1])
	print(f"Best RMS is {runs[0][1]} for run {runs[0][0]}.")
	
	plot_bars([x[1] for x in runs], [c[0] for c in param_candidates], "Add Parameter", "Parameter", "RMS")
	plt.tight_layout()
	if save_fname == None:
		plt.show()
	elif type(save_fname) == str:
		plt.savefig(save_fname)
	
	plt.close()
	
	return(runs)

def ommit_parameter(par_dict, lin_df, param_candidates, spfit_path=None, save_fname=None):
	runs = []
	orig_par_dict = par_dict.copy()
	
	param_candidates.insert(0, "Initial")

	
	for i, param in enumerate(param_candidates):
		par_dict =  orig_par_dict.copy()
		for j, cparam in enumerate(par_dict["PARAMS"]):
			if cparam[0] == param:
				del par_dict["PARAMS"][j]
				break
		
		results = run_spfit_v(par_dict, lin_df, spfit_path)
		RMS = parse_fit_result(results["message"])["RMS"]
		runs.append([i, RMS, par_dict["PARAMS"].copy()])

	runs = sorted(runs, key=lambda x: x[1])
	print(f"Best RMS is {runs[0][1]} for run {runs[0][0]}.")
	
	plot_bars([x[1] for x in runs], param_candidates, "Neglect Parameter", "Parameter", "RMS")
	plt.tight_layout()
	if save_fname == None:
		plt.show()
	elif type(save_fname) == str:
		plt.savefig(save_fname)
	
	plt.close()
	
	return(runs)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="pyfit_pycat is a wrapper around spfit and spcat for additional functionality")
	
	subparsers = parser.add_subparsers()
	
	sp = subparsers.add_parser("residuals")
	sp.set_defaults(action="residuals")
	sp.add_argument("FILE(s)", nargs="+", help="Specify Cat and Lin File")
	sp.add_argument("--QUERY", help="Querystring to filter transitions")
	
	sp = subparsers.add_parser("mixing")
	sp.set_defaults(action="mixing")
	sp.add_argument("FILE", help="Specify Egy File")
	sp.add_argument("QUERY", help="Querystring to filter energy levels")
	sp.add_argument("SAVENAME", nargs="?", help="Path to output file")
	
	args = vars(parser.parse_args())
	
	if args.get("action") == "residuals":
		files = args.get("FILE(s)")
		if len(files) == 1:
			cat_df = cat_to_df(f"{files[0]}.cat")
			lin_df = lin_to_df(f"{files[0]}.lin")
		else:
			cat_df = cat_to_df(f"{files[0]}")
			lin_df = lin_to_df(f"{files[1]}")
		residuals(cat_df, lin_df, args.get("QUERY"))
	
	elif args.get("action") == "mixing":
		print(args.get("FILE"))
		egy_df = egy_to_df(args.get("FILE"))
		mixing_coefficient(egy_df, args.get("QUERY"), args.get("SAVENAME"))

	else:
		print("No known action was specified")

		## Here you can write more complex stuff, that is not supported with the CLI interface
		## Some examples are right beneath (currently commented out)
		
		# var_dict = parvar_to_dict(r"path/to/your/project/molecule.var")
		# par_dict = parvar_to_dict(r"path/to/your/project/molecule.par")
		# int_dict = int_to_dict(r"path/to/your/project/molecule.int")
		# lin_df = lin_to_df(r"path/to/your/project/molecule.lin")
		# cat_df = cat_to_df(r"path/to/your/project/molecule.cat")
		# egy_df = egy_to_df(r"path/to/your/project/molecule.egy")
		
		## Best Candidate to add to Fit
		# cands = [[140101, 0.0, 1e+37], [410101, 0.0, 1e+37]]
		# add_parameter(par_dict, lin_df, cands, r"SPFIT_SPCAT")
		
		# Best Candidate to neglect from Fit
		# cands = [320101, 230101]
		# ommit_parameter(par_dict, lin_df, cands, r"SPFIT_SPCAT")
		
		## Plot Energies
		# energy_levels(egy_df, [1], range(10))
		
		## Check Crossings
		# check_crossings(egy_df, [1], range(10))
		
		## Plot Mixing Coefficients
		# mixing_coefficient(egy_df, "qn4 == 1 and qn2 < 20 and qn1 < 20 and qn1==qn2+qn3")

		## Residuals
		# residuals(cat_df, lin_df)