from __future__ import division
import numpy as np

def driftSim(n1 = 1, N = 100, s = 0):
    count = 1
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


def multipleDriftSims(n1 = 1, N = 100, s = 0, reps = 5):
    list_reps = []
    for i in range(reps):
        list_reps.append(driftSim(n1 = n1, N = N, s = s))
    return list_reps
