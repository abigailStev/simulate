import numpy as np
from astropy.table import Table, Column
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from scipy import fftpack
from multiprocessing import Pool, Array, Process
from scipy.optimize import brentq, least_squares

font_prop = font_manager.FontProperties(size=18)

def random_walk(n_steps=100, step_size=1.0):
    """
    Makes an array of random walk steps.

    Parameters
    ----------
    n_steps : int
        Number of steps in the random walk. Also the size of the output array.

    step_size : float
        Size of each step in the random walk.

    Returns
    -------
    path : np.array of floats
        The path of the random walk, aka the value at each step in the walk.
    """
    r = np.random.RandomState()
    path = np.zeros(n_steps)
    path[0] = 0.0
    for i in range(n_steps - 1):
        if (r.rand() >= 0.5):
            path[i + 1] = path[i] + step_size
        else:
            path[i + 1] = path[i] - step_size
            #         if path[i+1] >= np.pi or path[i+1] <= -np.pi:
            #             print i+1
            #             return i+1
    return path


def find_nearest(array, value):
    """
    Thanks StackOverflow!

    Parameters
    ----------
    array : np.array of ints or floats
        1-D array of numbers to search through. Should already be sorted from
        low values to high values.

    value : int or float
        The value you want to find the closest to in the array.

    Returns
    -------
    array[idx] : int or float
        The array value that is closest to the input value.

    idx : int
        The index of the array of the closest value.
    """
    idx = np.searchsorted(array, value, side="left")
    if idx == len(array) or np.fabs(value - array[idx - 1]) < \
            np.fabs(value - array[idx]):
        return array[idx - 1], idx - 1
    else:
        return array[idx], idx


def phase_angle(complex_number):
    return np.arctan2(complex_number.imag, complex_number.real)


def vecrotate(theta, complex_number):
    print "Theta, 4th segment:", theta[3]
    print "Before rotation, 4th segment"
    print "Abs:", np.abs(complex_number[3])
    print "Angle:", phase_angle(complex_number[3])
    x = complex_number.real
    y = complex_number.imag
    xrot = x * np.cos(theta) - y * np.sin(theta)
    yrot = x * np.sin(theta) + y * np.cos(theta)
    rotated_complex_number = xrot + yrot*1j
    # print type(rotated_complex_number)
    # print type(rotated_complex_number[0])
    # print np.shape(rotated_complex_number)
    print "After rotation, 4th segment"
    print "Abs:", np.abs(rotated_complex_number[3])
    print "Angle:", phase_angle(rotated_complex_number[3])
    return rotated_complex_number


def phils_way(ft, n_seg, ifund=2, iharm=4,):
    print "Shape ft:", np.shape(ft)
    print "n_seg: ", n_seg

    phi1obs = phase_angle(ft[ifund,:])
    phi2obs = phase_angle(ft[iharm,:])

    # phicor = -1. * phi2obs / 2.
    phicor = phi1obs

    # ftr2harm, fti2harm = vecrotate(phicor, ft.real[iharm], ft.imag[iharm])
    ft2fund = vecrotate(phicor, ft[ifund,:])

    print "Shape ft2fund:", np.shape(ft2fund)

    ft2harm = ft[iharm,:]
    phi1obs2 = phase_angle(ft2fund)
    phi2obs2 = phase_angle(ft2harm)

    obspsi = (phi2obs2 - phi1obs2) / 2. + np.pi / 4
    print "Shape obspsi:", np.shape(obspsi)
    # print obspsi
    fig, ax1 = plt.subplots(1, 1, figsize=(9, 5))
    bins_h, edges_h, patches_h = ax1.hist(obspsi, bins=70, range=[-3.142, 3.142],
                                          normed=False)
    ax1.set_xlabel("Psi (radians)")
    ax1.set_ylabel("Occurrences")
    plt.show()
    # plt.close()
    # ipsi = int(obspsi / (2. * np.pi / float(npsibins)))
    # psihist[ipsi] = psihist[ipsi] + 1.
    # psihist = 1

    cshr = ft2fund.real * ft2harm.real + ft2fund.imag * ft2harm.imag
    cshi = ft2fund.real * ft2harm.imag - ft2fund.imag * ft2harm.real
    print "Shape cshr:", np.shape(cshr)
    print "Shape cshi:", np.shape(cshi)
    csharm_re = np.mean(cshr)
    csharm_im = np.mean(cshi)
    # print csharm_re
    # print csharm_im
    csharm = csharm_re + 1j*csharm_im

    powharm = np.mean(np.abs(ft2harm))
    powfund = np.mean(np.abs(ft2fund))
    # print powharm
    # print powfund
    cohpsi = (csharm_re ** 2 + csharm_im ** 2) / (powharm * powfund)
    # print cohpsi
    # cohpsi = np.abs(csharm) / (powharm * powfund)

    # psifinal = -1. * ((np.arctan2(csharm_im, csharm_re)) / 2. + np.pi / 4.)
    psifinal = -1. * ((phase_angle(csharm) / 2.) + np.pi / 4.)

    if psifinal < 0:
        psifinal = psifinal + np.pi

    # My coherence is super large, so subtracting it from 1 makes this negative.
    errpsi = np.sqrt((1. - cohpsi) / (2. * cohpsi * float(n_seg)))

    print "Psi = ", psifinal, " +/- ", errpsi / 2.
    # exit()
    return


def fit_for_d(true_psi, psi_m):
    delta = np.abs(psi_m - true_psi)
    d_m = np.where(delta >= np.pi/2., np.pi - delta, delta)
    chisq = np.sum(d_m ** 2)
    return chisq

def for_each_h_test():
    # if h_offset < 0:
    #     print h_offset + np.pi
    # else:
    #     print h_offset
    n_seg = 12000
    num_bins = 64
    meta_dict = {'freq': 4.0,
                 # frequency of fundamental, in Hz (harmonic assumed to be 2*freq)
                 'dt': 0.0078125,  # time step between time bins, in seconds
                 'n_bins': num_bins * n_seg,
                 # int number of time bins in one segment
                 # 'amp1_ci': 100.,  # amplitude of fundamental of CI, in cts/s
                 # 'amp2_ci': 66.,  # amplitude of harmonic of CI, in cts/s
                 # 'mean_ci': 1000.,  # mean count rate of CI, in cts/s
                 'amp1_ref': 100.,  # amplitude of fundamental of ref, in cts/s
                 'amp2_ref': 66.,  # amplitude of harmonic of ref, in cts/s
                 'mean_ref': 1000.}  # mean count rate of ref, in cts/s

    # exposure = meta_dict['n_bins'] * meta_dict['dt']
    # print exposure

    tiny_bins = np.arange(0, meta_dict['n_bins'],
                          0.1)  # 10 tiny bins per 1 actual bin, to make a smooth sine wave
    period = 1.0 / meta_dict['freq']  # Period of sine waves, in seconds
    bpp = period / meta_dict['dt']  # Number of bins per period of sine wave

    ## How quickly the random walk reaches np.pi on average will set the Q-value of the QPO
    phase_walk = random_walk(n_steps=meta_dict['n_bins'], step_size=np.pi / 16.)
    phase_walk_tiny = np.repeat(phase_walk,
                                10)  # Defining phase_walk over tiny_bins

    # ci_fund = meta_dict['amp1_ci'] * np.sin(
    #     2.0 * np.pi * tiny_bins / bpp + phase_walk_tiny)
    # ci_harm = meta_dict['amp2_ci'] * np.sin(
    #     4.0 * np.pi * tiny_bins / bpp + 2 * (phase_walk_tiny + h_offset))
    ref_fund = meta_dict['amp1_ref'] * np.sin(
        2.0 * np.pi * tiny_bins / bpp + phase_walk_tiny)
    ref_harm = meta_dict['amp2_ref'] * np.sin(
        4.0 * np.pi * tiny_bins / bpp + 2 * (phase_walk_tiny + h_offset))

    # smooth_signal_ci = ci_fund + ci_harm + meta_dict['mean_ci']
    smooth_signal_ref = ref_fund + ref_harm + meta_dict['mean_ref']
    # signal_ci = np.mean(np.array_split(smooth_signal_ci, meta_dict['n_bins']),
    #                     axis=1)
    signal_ref = np.mean(np.array_split(smooth_signal_ref, meta_dict['n_bins']),
                         axis=1)
    # signal_ci[signal_ci < 0] = 0
    signal_ref[signal_ref < 0] = 0

    # noisy_signal_ci = signal_ci
    noisy_signal_ref = signal_ref

    meta_dict['n_bins'] = meta_dict['n_bins'] / n_seg

    # lc_ci = np.reshape(noisy_signal_ci, (n_seg, num_bins)).T
    lc_ref = np.reshape(noisy_signal_ref, (n_seg, num_bins)).T

    ## Initializations
    fourier = Table()
    fourier['FREQUENCY'] = Column(
        fftpack.fftfreq(meta_dict['n_bins'], d=meta_dict['dt']))
    fourier['POWER_CI'] = Column(np.zeros(meta_dict['n_bins']),
                                 dtype=np.float64)
    fourier['POWER_REF'] = Column(np.zeros(meta_dict['n_bins']),
                                  dtype=np.float64)
    fourier['CROSS'] = Column(np.zeros((meta_dict['n_bins'], n_seg)),
                              dtype=np.complex128)
    fourier['CROSS_AVG'] = Column(np.zeros(meta_dict['n_bins']),
                                  dtype=np.complex128)

    ## Subtracting the mean off each value of 'rate'
    mean_ref = np.mean(lc_ref, axis=0)
    rate_sub_mean_ref = np.subtract(lc_ref, mean_ref)
    fft_data_ref = fftpack.fft(rate_sub_mean_ref, axis=0)

    ##############
    ## PHIL'S WAY
    ##############

    # phils_way(fft_data_ref, n_seg, ifund=2, iharm=4)
    # conj = fft_data_ref[2, :] * np.conj(fft_data_ref[4, :])
    # mean_freqcross = np.mean(conj)
    # angle = phase_angle(mean_freqcross) % np.pi
    # tricky_diffs = np.append(tricky_diffs, angle)

    # Using equation 3 from Ingram and van der Klis 2015
    phi_h = phase_angle(fft_data_ref[4, :])
    phi_f = phase_angle(fft_data_ref[2, :])
    psi_m = ((phi_h - 2. * phi_f) / 2. % np.pi)  # psi per segment m
    # but multiplying phi_f by 2 and then dividing afterward, since phi_f is
    # better-defined than phi_h (per meeting notes with Phil, 7 Nov)

    # print fit_for_d(0, psi_m)
    # print fit_for_d(np.pi/2, psi_m)
    # print fit_for_d(-np.pi/2, psi_m)
    # print fit_for_d(np.pi, psi_m)

    many_psis = np.arange(0, 3.142, 0.01 * 3.142)
    # print "starting"
    many_chisqs = [fit_for_d(x, psi_m) for x in many_psis]
    # print "stopping"
    # results = brentq(fit_for_d, 0, np.pi, args=(psi))
    min_index = np.argmin(many_chisqs)
    # print "Min chisq:", many_chisqs[min_index]
    # print "True psi:", many_psis[min_index]

    deltas = np.abs(psi_m - many_psis[min_index])
    d_m = np.where(deltas >= np.pi / 2., np.pi - deltas, deltas)

    # Equation 5 from Ingram and van der Klis 2015
    # delta = np.where(psi >= np.pi/2., np.pi - psi, psi)

    fig, ax1 = plt.subplots(1, 1, figsize=(9, 5))
    bins_h, edges_h, patches_h = ax1.hist(psi_m, bins=150,
                                          range=[0, 3.142],
                                          normed=False)
    plt.close()
    # print edges_h[np.argmax(bins_h)]
    # print edges_h[np.argmax(bins_h)] - np.pi/4

    # There's a factor of pi/4 that needed to be subtracted to get the
    # answer back. Not clear where this comes from, but Phil found it too in
    # figuring out his method. Perhaps with how the FT is done?
    diffs = np.append(diffs, edges_h[np.argmax(bins_h)])
    true_psis = np.append(true_psis, many_psis[min_index])

    if h_offset % 1 == 0:
        print "\t", h_offset



if __name__ == "__main__":
    # all_h_offsets = np.arange(-3, 3, 0.05)
    # all_h_offsets = np.arange(-3, 3, 1)
    all_h_offsets = np.asarray([0.0, 0.3])
    diffs = np.asarray([])
    true_psis = np.asarray([])

    p = multiprocessing.Pool(len(all_h_offsets))

    print "Done!"

    out_tab = np.column_stack((all_h_offsets, diffs))
    np.savetxt("psi_maxhistogram.txt", out_tab)

    out_tab = np.column_stack((all_h_offsets, true_psis))
    np.savetxt("fitting_for_true_psis.txt", out_tab)


    fig, ax = plt.subplots(1,1,figsize=(7,7))
    ax.scatter(all_h_offsets % np.pi, diffs)
    ax.set_xlim(-0.1, 3.1)
    ax.set_ylim(-0.1, 3.1)
    ax.grid(b=True, which='major', color='gray', linestyle='-')
    ax.set_title("Bin of max histogram value")
    ax.set_xlabel("Original harmonic offset, mod pi", fontproperties=font_prop)
    ax.set_ylabel("Measured phase difference, mod pi", fontproperties=font_prop)
    plt.savefig("psi_maxhistogram.png", dpi=300)
    plt.close()
    print "psi_maxhistogram.png"

    fig, ax = plt.subplots(1,1,figsize=(7,7))
    ax.scatter(all_h_offsets % np.pi, true_psis)
    ax.set_xlim(-0.1, 3.1)
    ax.set_ylim(-0.1, 3.1)
    ax.grid(b=True, which='major', color='gray', linestyle='-')
    ax.set_title("Fitting for true psi")
    ax.set_xlabel("Original harmonic offset, mod pi", fontproperties=font_prop)
    ax.set_ylabel("Measured phase difference, mod pi", fontproperties=font_prop)
    plt.savefig("fitting_for_true_psis.png", dpi=300)
    plt.close()
    print "fitting_for_true_psis.png"