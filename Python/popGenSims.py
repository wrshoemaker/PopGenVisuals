from __future__ import division
import numpy as np

def driftSim(n1 = 1, N = 100, s = 0):
    count = 0
    gen_list = []
    freq_list = []
    while 0 < n1 < N:
        n0 = N - n1
        p1 = n1 * (1 + s) / (n0 + (n1 * (1 + s)))
        n1 = np.random.binomial(N, p1)
        count += 1
        gen_list.append(count)
        freq_list.append(p1)
    return (gen_list, freq_list)


def driftSimBottle(n1 = 1, N = 10000, s = 0, bottleReduct = 100):
    count = 1
    gen_list = [1]
    freq_list = [n1/N]
    while 0 < count < 80:
        n0 = N - n1
        p1 = n1 * (1 + s) / (n0 + (n1 * (1 + s)))
        n1 = np.random.binomial(N, p1)
        count += 1
        gen_list.append(count)
        freq_list.append(p1)
    n1_bottle = round(n1/bottleReduct)
    N_bottle = round(N/bottleReduct)
    while 0 < n1_bottle < N_bottle:
        n0_bottle = N_bottle - n1_bottle
        p1_bottle = n1_bottle * (1 + s) / (n0_bottle + (n1_bottle * (1 + s)))
        n1_bottle = np.random.binomial(N_bottle, p1_bottle)
        count += 1
        gen_list.append(count)
        freq_list.append(p1_bottle)
    return (gen_list, freq_list)

def multipleDriftBottleSims(n1 = 1, N = 100, s = 0, bottleReduct = 100, reps = 5):
    list_reps = []
    count = reps
    while count > 0:
        out = driftSimBottle(n1 = n1, N = N, s = s, bottleReduct = bottleReduct)
        if out[0][-1] > 80:
            list_reps.append(out)
            count -= 1
    return list_reps

def multipleDriftSims(n1 = 1, N = 100, s = 0, reps = 5):
    list_reps = []
    count = 5
    for i in range(reps):
        list_reps.append(driftSim(n1 = n1, N = N, s = s))
    return list_reps
