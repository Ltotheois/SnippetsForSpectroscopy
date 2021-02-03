# FFT-Correction

*fft_correction.py* is a script to eliminate the baseline from measured spectra by using fast Fourier transform. This file is also available as an *.exe* file under releases.

On top a slider can be used to set the cut-off frequency for the FFT-coefficients, which can be seen in the first plot. Also clicking in the coefficient plot will change the cut-off frequency. The cut-off frequency is indicated by a dashed line in the FFT-coefficients plot. If keyboard use is preffered, the arrow keys (in combination with shift) can be used to change the cut-off frequency.

Beneath the original spectrum and the current calculated baseline are shown, with the resulting spectrum beneath.

Files are selected in the initial pop-up window. If a matching *.info* file is given, the duration of the fit is calculated correctly, otherwise it is set to 30 seconds. The duration is only important for the x-axis of the FFT-coefficient plot.

Press enter or space to save the changes for a file and go to the next file, press esc to skip the file, press ctrl+q to close the app. Press ctrl+left to go to the previous file without saving.