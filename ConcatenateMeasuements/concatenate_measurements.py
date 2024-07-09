#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Luis Bonah
# Description: Short script to concatenate measurements


import glob
import numpy as np
import argparse


def concat_measurements(filenames, outputfilename):
	if not filenames:
		print('No files were selected.')
		return

	datasets = [np.genfromtxt(filename, delimiter='\t') for filename in filenames]
	shapes = np.array([dataset.shape for dataset in datasets])
	
	if len(shapes.shape) != 2:
		raise ValueError('The shape of the data was not understood.')

	if not np.all(shapes[0, 1] == shapes[:, 1]):
		raise ValueError('Files have different number of columns.')
	
	xmins = np.array([np.nanmin(dataset[:, 0]) for dataset in datasets])
	xmaxs = np.array([np.nanmax(dataset[:, 0]) for dataset in datasets])

	sort_mask_mins = np.argsort(xmins)
	sort_mask_maxs = np.argsort(xmaxs)

	sort_masks_are_equivalent = (sort_mask_mins == sort_mask_maxs).all()

	if not sort_masks_are_equivalent:
		raise ValueError('')
	
	xmins = xmins[sort_mask_mins]
	xmaxs = xmaxs[sort_mask_maxs]

	if not (xmaxs[:-1] <= xmins[1:]).all():
		print('The measurements are overlapping. Data will be interwoven. Please check, if this is desired.')
		complete_data = np.concatenate(datasets)
		sort_mask = np.argsort(complete_data[:, 0])
		complete_data = complete_data[sort_mask]
	else:
		complete_data = np.concatenate([datasets[i] for i in sort_mask_mins])
	
	if not outputfilename:
		outputfilename = input('Please specify the output file: ')
	
	np.savetxt(outputfilename, complete_data, delimiter='\t')
	print(f'The concatenated spectrum was written to \'{outputfilename}\'.')



if __name__ == '__main__':
	epilog = 'Concatenate multiple measurements into a single file'
	parser = argparse.ArgumentParser(prog='Concatenate Measurements', epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
	
	parser.add_argument('files', nargs='*', type=str, help='Glob string for files to be loaded')
	parser.add_argument('--output', '-o', type=str, default=None, help='Output file')

	args = parser.parse_args()
	
	
	if args.files:
		files = [file for glob_str in args.files for file in glob.glob(glob_str, recursive=True)]
	else:
		files = []
		
	outputfilename = args.output
	concat_measurements(files, outputfilename)
