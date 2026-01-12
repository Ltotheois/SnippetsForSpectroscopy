#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import pyckett
import requests
import uncertainties
import argparse
import pandas as pd
import numpy as np
from datetime import datetime
from scipy import optimize, special
import matplotlib
import matplotlib.pyplot as plt
import shutil

matplotlib.rcParams['axes.formatter.useoffset'] = False


BOLD = "\033[1m"
RESET = "\033[0m"
GREEN = "\033[92m"
RED = "\033[91m"


molecule_identifiers = {
    'OCS': 60503,
    'OCS, v2=1': 60504,
    'O13CS': 61502,
    'OC33S': 61503,
    '17OCS': 61504,
    'OC34S': 62505,
    '18OCS': 62506,
}

molecule_labels = {label: key for key, label in molecule_identifiers.items()}

def get_cat_file_from_cdms(label, url, timeout=30):
    """
    Downloads a cat file from a URL and returns a dataframe.
    Raises requests.exceptions.RequestException on error.
    """
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    data = response.content.decode("utf-8", errors="ignore")
    cat = pyckett.cat_to_df(pyckett.str_to_stream(data))
    cat['filename'] = label
    return cat

def Voigt(derivative, x, x0, amp, fwhm_gauss, fwhm_lorentz, *poly_vals):
    sigma = fwhm_gauss / (2 * np.sqrt(2 * np.log(2)))
    gamma = fwhm_lorentz / 2
    if gamma == sigma == 0:
        return [0 if i != x0 else np.inf for i in x]

    z = (x - x0 + 1j * gamma) / (sigma * np.sqrt(2))
    wz = special.wofz(z)
    w0 = special.wofz((1j * gamma) / (sigma * np.sqrt(2)))

    if derivative == 0:
        tmp = lambda x, x0, wz, sigma, gamma: np.real(wz) / (sigma * np.sqrt(2 * np.pi))
    elif derivative == 1:
        tmp = (
            lambda x, x0, wz, sigma, gamma: 1
            / (sigma**3 * np.sqrt(2 * np.pi))
            * (gamma * np.imag(wz) - (x - x0) * np.real(wz))
        )
    elif derivative == 2:
        tmp = (
            lambda x, x0, wz, sigma, gamma: 1
            / (sigma**5 * np.sqrt(2 * np.pi))
            * (
                gamma * (2 * (x - x0) * np.imag(wz) - sigma * np.sqrt(2 / np.pi))
                + (gamma**2 + sigma**2 - (x - x0) ** 2) * np.real(wz)
            )
        )
    else:
        raise NotImplementedError(
            "Only the zeroth, first, and second derivatives of a Gaussian are implemented."
        )

    ys = tmp(x, x0, wz, sigma, gamma)
    if derivative == 1:
        xspan = max(sigma, gamma) * 2
        tmp_xs = np.linspace(-xspan, +xspan, 1000)
        ymax = np.max(tmp(tmp_xs, 0, w0, sigma, gamma))
    else:
        ymax = tmp(0, 0, w0, sigma, gamma)
    if ymax == 0:
        ymax = 1
    ys *= amp / ymax
    ys += np.polyval(poly_vals, x)

    return ys

def fit_lineshape(lineshape, xs, ys, baselinerank=0):
    xmin, xmax = xs.min(), xs.max()
    x0 = (xmin + xmax) / 2
    yptp = np.ptp(ys)
    y0 = 0

    w0 = (xmax - xmin) / 10
    wmin = 0
    wmax = (xmax - xmin)

    amp_min, amp_max = -3 * yptp, 3 * yptp

    p0 = [x0, y0, w0, w0] + [0] * baselinerank
    bounds = [
        [xmin, amp_min, wmin, wmin] + [-np.inf] * baselinerank, 
        [xmax, amp_max, wmax, wmax] + [+np.inf] * baselinerank, 
    ]

    popt, pcov = optimize.curve_fit(lineshape, xs, ys, p0=p0, bounds=bounds)
    perr = np.sqrt(np.diag(pcov))
    return(popt, perr)

def run_calibration(cat_df, xmin=None, xmax=None, max_x_deviation=0.1, existing_measurements=None, skip_figure=False, folder_path=None, individual_figure=False, baselinerank=0):
    now = datetime.now()
    timestamp = now.strftime("%Y/%m/%d %H:%M:%S")
    timestamp_folder = now.strftime("%Y%m%d_%H%M%S")
    if not folder_path:
        folder_path = 'Calibration_' + timestamp_folder
    os.makedirs(folder_path, exist_ok=True)

    listfname = os.path.join(folder_path, 'Calibration.list')
    measurements_folder = os.path.join(folder_path, 'Measurements')
    os.makedirs(measurements_folder, exist_ok=True)
    
    query = []
    if xmin:
        query.append(f' (x > {xmin})')
    if xmax:
        query.append(f' (x < {xmax})')
    
    if query:
        query = ' and '.join(query)
        cat_df = cat_df.query(query)

    if not len(cat_df):
        print('⚠️  No calibration lines found in your frequency range.')
        return
    
    list_content = ['probe.center'] + [f'{x:12.4f}' for x in cat_df['x']]
    list_content = '\n'.join(list_content)

    with open(listfname, 'w+') as file:
        file.write(list_content)


    print(f'✅ Wrote {len(cat_df)} transitions to the file \'{listfname}\'.\n')
    print(f'{BOLD}1{RESET} Run these measurements')
    print(f'{BOLD}2{RESET} Copy the measurements into the folder \'{measurements_folder}\'')
    resp = input(f'{BOLD}3{RESET} Choose if you want to analyze the measurements (Y/n)')

    if resp.lower() == 'n':
        return
    
    if existing_measurements:
        for file in glob.glob(existing_measurements):
            shutil.copy(file, measurements_folder)

    glob_string = os.path.join(measurements_folder, '*.dat')
    measurement_fnames = glob.glob(glob_string)

    while not len(measurement_fnames):
        resp = input('⚠️  No measurement files found, do you want to retry? (Y/n)')
        if resp.lower() == 'n':
            return
        
        glob_string = os.path.join(measurements_folder, '*.dat')
        measurement_fnames = glob.glob(glob_string)

    lineshape = lambda *args: Voigt(2, *args)
    fit_results = []
    labels = []

    for measurement_fname in measurement_fnames:
        data = np.genfromtxt(measurement_fname, delimiter='\t')
        xs, ys = data[:, 0], data[:, 1]
        exp_range= xs.min(), xs.max()
        xoffset = sum(exp_range)/2

        popt, perr = fit_lineshape(lineshape, xs-xoffset, ys, baselinerank=baselinerank)

        x0 = popt[0] + xoffset

        closest_idx = (cat_df["x"] - x0).abs().idxmin()
        closest_row = cat_df.loc[closest_idx]

        if abs(x0 - closest_row['x']) > max_x_deviation:
            print(f'⚠️  File was skipped due to large deviation from literature position: {measurement_fname}')
            continue

        x0_exp = uncertainties.ufloat(x0, perr[0])
        x0_lit = uncertainties.ufloat(closest_row['x'], closest_row['error'])

        fit_results.append((x0_exp.n, x0_lit.n, x0_exp.s, x0_lit.s))
        labels.append(closest_row['filename'])

        if individual_figure:
            fit_xs = np.linspace(*exp_range, 1000)
            fit_ys = lineshape(fit_xs-xoffset, *popt)

            fig, ax = plt.subplots()
            ax.plot(xs, ys, color="#DC267F", label='Experiment', linewidth=2)
            ax.plot(fit_xs, fit_ys, color="#648FFF", label='Fit', linewidth=1)

            ax.set_xlabel('Frequency [MHz]')
            ax.set_ylabel('Intensity [A.U.]')
            ax.set_xlim(*exp_range)

            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5))
            ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5))
            ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
            ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
            ax.grid(which='both', linestyle=':', linewidth=0.5, color='gray')

            ax.legend(loc='upper right')
            plt.tight_layout()
            plt.savefig(os.path.splitext(measurement_fname)[0] + ".png", dpi=600)
            plt.close()


    if not len(fit_results):
        print(f'⚠️  No file was analyzed!')
        return
    
    fit_results = np.array(fit_results)
    xs_exp = fit_results[:, 0]
    xs_lit = fit_results[:, 1] 
    us_exp = fit_results[:, 2]
    us_lit = fit_results[:, 3]
    dxs = xs_exp - xs_lit
    udxs = np.sqrt(us_exp**2 + us_lit**2)

    print('\n\n')
    report_string = [
        f'{BOLD}CALIBRATION REPORT{RESET}\n',
        f'Summary for calibration from {timestamp}\n',
    ]

    report_string.append('|   Exp [MHz]   |   Lit [MHz]   | ΔExp [MHz] | ΔLit [MHz] |  Dev [kHz] | ΔDev [kHz] |           Label           |')
    report_string.append('|---------------|---------------|------------|------------|------------|------------|---------------------------|')
    for label, xexp, xlit, uexp, ulit in zip(labels, xs_exp, xs_lit, us_exp, us_lit):
        diff_unc = np.sqrt(uexp**2 + ulit**2) * 1000
        diff = (xexp-xlit)*1000

        if diff_unc < abs(diff):
            color = RED
            color_reset = RESET
        else:
            color = GREEN
            color_reset = RESET

        label = label[:25]
        report_string.append(f'| {xexp:13.4f} | {xlit:13.4f} | {uexp:10.4f} | {ulit:10.4f} | {color}{diff:10.2f}{color_reset} | {color}{diff_unc:10.2f}{color_reset} | {label:25} |')

    report_string.append(f'\nUsed a total of {len(xs_exp)} transitions')
    report_string.append(f'\nDeviations RMS: {np.sqrt(np.mean(dxs**2))*1000:10.2f} kHz')

    i_max = np.argmax(np.abs(dxs))
    tmp = dxs[i_max]
    report_string.append(f'Max Deviation: {tmp*1000:11.2f} kHz (for x = {xs_exp[i_max]:.2f})')

    i_rel_max = np.argmax(np.abs(dxs/udxs))
    tmp = dxs[i_rel_max] / udxs[i_rel_max]
    report_string.append(f'Max Rel Deviation: {tmp:7.2f}     (for x = {xs_exp[i_max]:.2f})\n')


    report_string = '\n'.join(report_string)
    print(report_string)
    with open(os.path.join(folder_path, 'Report.txt'), 'w+', encoding='utf-8') as file:
        ansi_escape = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')
        clean_text = ansi_escape.sub('', report_string)
        file.write(clean_text)

    fig, ax = plt.subplots()
    ax.errorbar(xs_exp, dxs, udxs, marker='o', linestyle='none', color='#648FFF')
    ax.axhline(y=0, lw=0.5, color='gray')

    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Residuals [MHz]')

    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5))
    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
    ax.grid(which='both', linestyle=':', linewidth=0.5, color='gray')

    plt.tight_layout()
    plt.savefig(os.path.join(folder_path, 'Figure.png'), dpi=600)

    if not skip_figure:
        plt.show()

    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calibrate frequency accuracy with literature data"
    )

    parser.add_argument(
        "molecules",
        type=str,
        nargs="*",
        help="CDMS molecule identifiers (e.g., 60503) or path to *.cat file"
    )

    parser.add_argument(
        '-r',
        "--range",
        type=float,
        nargs=2,
        help="Specify frequency range",
        default=(None, None),
    )

    parser.add_argument(
        '-m',
        '--measurements',
        type=str,
        help='Glob string to already existing measurements',
    )

    parser.add_argument(
        '-s',
        '--skip',
        action="store_true",
        help="Do not automatically open the summary figure"
    )

    parser.add_argument(
        '-f',
        '--folderpath',
        type=str,
        default = None,
        help="Folder to save calibration to"
    )

    parser.add_argument(
        '-i',
        '--individual',
        action="store_true",
        help="Create figures for individual lines"
    )

    parser.add_argument(
        '-b',
        '--baselinerank',
        type=int,
        default=0,
        help="Rank of baseline polynom for fit of individual lines"
    )

    args = parser.parse_args()

    cat_dfs = []
    if not len(args.molecules):
        molecule_urls = {key: f"https://cdms.astro.uni-koeln.de/classic/entries/c{id:06.0f}.cat" for key, id in molecule_identifiers.items()}
        cat_dfs = [get_cat_file_from_cdms(label, url) for label, url in molecule_urls.items()]
    else:
        molecule_urls = {}
        for id in args.molecules:
            if os.path.isfile(id):
                tmp = pyckett.cat_to_df(id)
                tmp['filename'] = os.path.basename(id)
                cat_dfs.append(tmp)
            else:
                id = int(id)
                label = molecule_labels.get(id, id)
                url = f"https://cdms.astro.uni-koeln.de/classic/entries/c{id:06.0f}.cat"
                cat_dfs.append(get_cat_file_from_cdms(label, url))
    cat_df = pd.concat(cat_dfs).reset_index(drop=True)

    run_calibration(cat_df, *args.range, existing_measurements=args.measurements, folder_path=args.folderpath, skip_figure=args.skip, individual_figure=args.individual, baselinerank=args.baselinerank)
