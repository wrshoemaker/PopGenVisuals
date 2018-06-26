from __future__ import division
import os, math
import numpy as np
import  matplotlib.pyplot as plt

mydir = os.path.expanduser("~/GitHub/PopGenVisuals/")

def make_plot(N):
    # s_N is a value, s_S will be a range
    #(s_N / s_S)
    s_S = np.logspace(-6, -1, num = 1000, endpoint=True, base=10.0)
    numer_S = [(1 - np.exp(-2*N*x)) for x in s_S]
    s_div_s_1 = 0.01 / s_S
    s_div_s_2 = 0.001 / s_S
    s_div_s_3 = 0.0001 / s_S

    denom_N_1 = 1 - np.exp( -2*N*0.01)
    denom_N_2 = 1 - np.exp( -2*N*0.001)
    denom_N_3 = 1 - np.exp( -2*N*0.0001)
    y_1 = (numer_S / denom_N_1) * s_div_s_1
    y_2 = (numer_S / denom_N_2) * s_div_s_2
    y_3 = (numer_S / denom_N_3) * s_div_s_3

    fig = plt.figure()
    plt.plot(s_S, y_1, lw = 2, c = '#FF6347', label = r'$s_{N}=0.01$')
    plt.plot(s_S, y_2, lw = 2, c = '#FFA500', label = r'$s_{N}=0.001$')
    plt.plot(s_S, y_3, lw = 2, c = '#87CEEB', label = r'$s_{N}=0.0001$')
    #plt.ylim(0, 1)
    plt.legend(loc = 'upper right')
    plt.xlabel(r'$s_{S}$', fontsize=24)
    plt.ylabel(r'$d_{N}/d_{S}$' , fontsize=16)
    plt.xscale('log')
    plt.yscale('log')
    fig.savefig(mydir + 'Figures/dn_ds.png', bbox_inches='tight', pad_inches = 0.4, dpi = 600)
    plt.close()

make_plot(10000)
