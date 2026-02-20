#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Luis Bonah
# Description: Short script to remove lines from measurements

import numpy as np
import matplotlib.pyplot as plt
import pyckett
import argparse


def remove_lines_from_spectrum(exp_xs, exp_ys, remove_xs, save_fname, linewidth=5, zerovalue=np.nan):
	exp_xs, exp_ys = exp.T
	exp_min, exp_max = exp_xs.min(), exp_xs.max()
	exp_ys_mod = exp_ys.copy()

	start_indices = np.searchsorted(exp_xs, remove_xs - linewidth, side="left")
	stop_indices = np.searchsorted(exp_xs, remove_xs + linewidth, side="right")

	for i_start, i_stop in zip(start_indices, stop_indices):
		exp_ys_mod[i_start:i_stop] = zerovalue

	data_out = np.vstack((exp_xs, exp_ys_mod)).T
	np.savetxt(save_fname, data_out, delimiter="\t")

if __name__ == '__main__':
	epilog = 'Remove known lines from an experimental spectrum'
	parser = argparse.ArgumentParser(prog='Remove Lines', epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
	
	parser.add_argument('exp_fname', type=str, help='Path to experimental file')
	parser.add_argument('lin_fname', type=str, help='Path to *.lin file with known frequencies')
	parser.add_argument('--output', '-o', type=str, default='SpectrumModified.csv', help='Output file')
	parser.add_argument('--linewidth', '-l', type=float, default=5, help='Distance to remove spectrum around line centers')

	args = parser.parse_args()
	
	linewidth = args.linewidth
	zerovalue = np.nan

	save_fname = args.output
	exp_fname = args.exp_fname
	lin_fname = args.lin_fname

	exp = np.genfromtxt(exp_fname)
	exp_xs, exp_ys = exp[:, 0], exp[:, 1]
	
	lin = pyckett.lin_to_df(lin_fname)
	remove_xs = lin["x"].values
	
	remove_lines_from_spectrum(exp_xs, exp_ys, remove_xs, save_fname, linewidth=linewidth, zerovalue=zerovalue)