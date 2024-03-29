# FFT-Correction

*fft_correction.py* is a script to eliminate the baseline from measured spectra by using a fast Fourier transform for real input.

Load files by selecting them in the initial pop-up window or by pressing *Replace Files* (ctrl+o).

The cut-off frequency for the FFT-coefficients can be set by moving the slider on top or clicking in the FFT-coefficientsplot. The cut-off frequency is indicated by a dashed line in the FFT-coefficients plot. If keyboard use is preferred, the arrow keys (additionally holding shift increases the step size) can be used to change the cut-off frequency.

Beneath are the original spectrum and the current calculated baseline shown, with the resulting spectrum on the bottom.

Press *Save* (enter or space) to save the FFT-corrected version of a file and go to the next file. The file will be saved with an additional "FFT" suffix before the file extension. Therefore do not keep files with such a name in the same folder, as they might be overwritten (as an example: correcting "1.dat" will overwrite (if exisiting) "1FFT.dat").

Press *Next* (Esc or ctrl+right) to go to the next file without saving, press *Previous* (ctrl+left) to go to the previous file without saving. Press *Quit* (ctrl+q) to close the script. 

You can zoom in the plots by selecting a range with the right mouse button. To go back to the initial zoom press *Reset Zoom* (ctrl+r). To save the currently displayed figure press *Save Figure* (ctrl+s).

*NEW*: With HPR (standing for high-pass filter rank) you can change the rank of the high-pass filter. Setting it to 0 just zeroes the fft-bins beneath the cutoff frequency.