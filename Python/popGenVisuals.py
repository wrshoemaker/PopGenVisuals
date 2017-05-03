from __future__ import division
import os, math
import popGenSims as pgs
import numpy as np
import  matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import animation
import pandas as pd
from scipy import stats

mydir = os.path.expanduser("~/GitHub/PopGenVisuals/")


def driftFig(N = 1000, n1 = 50,  s = 0, reps = 50, x_range = 500):
    sims = pgs.multipleDriftSims(n1 = n1, N = N, s = s, reps = reps)
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(1, 1, 1)
    plt.xlabel('Generation', fontsize=20)
    if s != 0:
        plt.ylabel('Frequency of beneficial allele', fontsize=20)
    else:
        plt.ylabel('Allele frequency', fontsize=20)
    for sim in sims:
        plt.plot(sim[0], sim[1], c = 'b', alpha = 0.5)
    fig.tight_layout()
    ax.set_xlim(0, x_range)
    fig_name = mydir + 'Figures/Drift/Drift_N_' + str(N) + '_n_'+ \
        str(n1) + '_s_' + str(s) + '.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def timeToLoss(p = 0.5):
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(1, 1, 1)
    plt.xlabel('Population size', fontsize=20)
    plt.ylabel('Time to loss', fontsize=20)
    Ns = np.logspace(1, 4, num = 1000, base=10.0)
    T_fix = [ ((-4 * N * p) / (1-p)) * math.log(p) for N in Ns]
    plt.plot(Ns, T_fix, c = 'b')
    fig.tight_layout()
    ax.set_xlim(0, 350)
    fig_name = mydir + 'Figures/Drift/time_to_loss.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def lossHetero(h_0 = 0.5, t = 500):
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(1, 1, 1)
    plt.xlabel('Generation', fontsize=20)
    plt.ylabel('Heterozygosity', fontsize=20)
    t_range = range(t)
    colors = ['#FF6347', '#FFA500', '#87CEEB']
    h_10 = [h_0 * math.exp(-t_i / (2 * 10)) for t_i in t_range]
    h_100 = [h_0 * math.exp(-t_i / (2 * 100)) for t_i in t_range]
    h_1000 = [h_0 * math.exp(-t_i / (2 * 1000)) for t_i in t_range]
    plt.plot(t_range, h_10, c = '#FF6347', lw = 2, label = 'N = 10')
    plt.plot(t_range, h_100, c = '#FFA500', lw = 2, label = 'N = 100')
    plt.plot(t_range, h_1000, c = '#87CEEB', lw = 2, label = 'N = 1000')
    ax.legend(loc='upper left', prop={'size':12})
    fig.tight_layout()
    fig_name = mydir + 'Figures/Drift/loss_hetero.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def fixProb(Nmax = 1000000, s = 0.01):
    def kimura(N, s):
        return (np.expm1(-2 * s) ) / (np.expm1(-2 * s * N))
    Ns = np.logspace(1, np.log10(Nmax), num = 1000, base=10.0)
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(1, 1, 1)
    plt.xlabel('Population size (N)', fontsize=20)
    plt.ylabel('Probability of fixation', fontsize=20)
    Pfix = [kimura(x, s) for x in Ns]
    plt.plot(Ns, Pfix, c = 'b', lw = 2)
    fig.tight_layout()
    #ax.set_xlim(0, 350)
    ax.set_xscale('log')
    fig_name = mydir + 'Figures/Drift/prob_fixation.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def galtonHist():
    IN = pd.read_csv(mydir + 'data/Galton.csv', sep = ',')
    M = IN.loc[IN['Gender'] == 'M'].Height
    F = IN.loc[IN['Gender'] == 'F'].Height
    fig = plt.figure()
    plt.hist(M, 20, fc='#87CEEB', histtype='bar', alpha=0.5, normed = False, label = 'Men')
    plt.hist(F, 20, fc='#FF6347', histtype='bar', alpha=0.5, normed = False, label = 'Women')
    plt.xlabel('Height in inches', fontsize = 18)
    plt.ylabel('Number of individuals', fontsize = 18)
    plt.ylim(0, 75)
    plt.legend(loc = 'upper left')
    fig.tight_layout()
    fig_name = mydir + 'Figures/galtonHist.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def galtonRegress():
    IN = pd.read_csv(mydir + 'data/Galton.csv', sep = ',')
    IN['Midparent'] = IN[['Father', 'Mother']].mean(axis=1)
    x_M = IN.loc[IN['Gender'] == 'M'].Midparent
    x_F = IN.loc[IN['Gender'] == 'F'].Midparent
    y_M = IN.loc[IN['Gender'] == 'M'].Height
    y_F = IN.loc[IN['Gender'] == 'F'].Height
    fig = plt.figure()
    plt.scatter(x_M, y_M, c='#87CEEB', marker='o', label='Men')
    slope_M, intercept_M, r_value_M, p_value_M, std_err_M = stats.linregress(x_M,y_M)
    predict_y_M = intercept_M + slope_M * x_M
    pred_error_M = y_M - predict_y_M
    degrees_of_freedom_M = len(x_M) - 2
    residual_std_error_M = np.sqrt(np.sum(pred_error_M**2) / degrees_of_freedom_M)
    plt.plot(x_M, predict_y_M, 'k-', lw = 5, c = 'black', label = '_nolegend_' )
    plt.plot(x_M, predict_y_M, 'k-', lw = 2, c = '#87CEEB', label = '_nolegend_')

    plt.scatter(x_F, y_F, c='#FF6347', marker='o', label='Women')
    slope_F, intercept_F, r_value_F, p_value_F, std_err_F = stats.linregress(x_F,y_F)
    predict_y_F = intercept_F + slope_F * x_F
    pred_error_F = y_F - predict_y_F
    degrees_of_freedom_F = len(x_F) - 2
    residual_std_error_F = np.sqrt(np.sum(pred_error_F**2) / degrees_of_freedom_F)
    plt.plot(x_F, predict_y_F, 'k-', lw = 5, c = 'black', label = '_nolegend_' )
    plt.plot(x_F, predict_y_F, 'k-', lw = 2, c = '#FF6347', label = '_nolegend_')

    plt.xlabel('Mid-parent height (inches)', fontsize = 18)
    plt.ylabel('Offspring height (inches)', fontsize = 18)
    plt.legend(loc = 'upper left')
    fig.tight_layout()
    fig_name = mydir + 'Figures/galtonRegress.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def multiLocusTraits(L = 2):
    freqs = [0.25, 0.5, 0.25]
    fig = plt.figure()
    width = .5
    ind = np.arange(len(freqs))
    h = plt.bar(ind, freqs, width=width, align='center', color = '#87CEEB', alpha = 0.8)
    plt.xticks(ind , freqs)
    fig.text(0.53, 0.45, 'Aa',  fontsize=18,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.24, 0.45, 'aa',  fontsize=18,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.83, 0.45, 'AA',  fontsize=18,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    xticks_pos = [0.65*patch.get_width() + patch.get_xy()[0] for patch in h]
    plt.xticks(xticks_pos, freqs,  ha='right', rotation=45)
    plt.xticks([0,1,2], [0,1,2], rotation='horizontal')
    plt.xlabel('Number of height increasing alleles', fontsize = 18)
    plt.ylabel('Phenotype frequency', fontsize = 18)
    fig.tight_layout()
    fig_name = mydir + 'Figures/multiLocusTraits/L1.png'
    fig.savefig(fig_name, bbox_inches = None, pad_inches = 0.8, dpi = 600)
    plt.close()


    freqs = [1/9, 2/9, 3/9, 2/9, 1/9]
    fig = plt.figure()
    width = 0.65
    ind = np.arange(len(freqs))
    h = plt.bar(ind, freqs, width=width, align='center', color = '#87CEEB', alpha = 0.8)
    plt.xticks(ind , freqs)
    xticks_pos = [0.65*patch.get_width() + patch.get_xy()[0] for patch in h]
    fig.text(0.26 , 0.35, 'aabb',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')

    fig.text(0.40, 0.40, 'Aabb',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.40, 0.35, 'aaBb',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')

    fig.text(0.545, 0.45, 'AAbb',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.545, 0.40, 'AaBb',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.545, 0.35, 'aaBB',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')

    fig.text(0.685, 0.40, 'AABb',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.685, 0.35, 'AaBB',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')

    fig.text(0.83 , 0.35, 'AABB',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')


    plt.xticks(xticks_pos, freqs,  ha='right', rotation=45)
    plt.xticks([0,1,2, 3, 4], [0,1,2, 3, 4], rotation='horizontal')
    plt.xlabel('Number of height increasing alleles', fontsize = 18)
    plt.ylabel('Phenotype frequency', fontsize = 18)
    fig.tight_layout()
    fig_name = mydir + 'Figures/multiLocusTraits/L2.png'
    fig.savefig(fig_name, bbox_inches = None, pad_inches = 0.8, dpi = 600)
    plt.close()


    freqs = [1/27, 3/27, 6/27,  7/27, 6/27, 3/27, 1/27]
    fig = plt.figure()
    width = 0.84
    ind = np.arange(len(freqs))
    h = plt.bar(ind, freqs, width=width, align='center', color = '#87CEEB', alpha = 0.8)
    plt.xticks(ind , freqs)
    xticks_pos = [0.65*patch.get_width() + patch.get_xy()[0] for patch in h]

    fig.text(0.22 , 0.19, 'aabbcc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')

    fig.text(0.33, 0.19, 'Aabbcc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.33, 0.24, 'aaBbcc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.33, 0.29 , 'aabbCc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')

    fig.text(0.4365, 0.19, 'AAbbcc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.4365, 0.24, 'AaBbcc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.4365, 0.29, 'AabbCc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.4365, 0.34, 'aaBBcc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.4365, 0.39, 'aaBbCc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.4365, 0.44, 'aabbCC',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')

    fig.text(0.544, 0.19, 'AABbcc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.544, 0.24, 'AaBbcc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.544, 0.29, 'AabbCc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.544, 0.34, 'aaBBcc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.544, 0.39, 'aaBbCc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.544, 0.44, 'aabbCC',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.544, 0.49, 'aabbCC',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')

    fig.text(0.655, 0.19, 'AABBcc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.655, 0.24, 'AABbCc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.655, 0.29, 'AAbbCC',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.655, 0.34, 'AaBBCc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.655, 0.39, 'AaBbCC',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.655, 0.44, 'aaBBCC',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')

    fig.text(0.762, 0.19, 'AABBCc',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.763, 0.24, 'AABbCC',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.764, 0.29 , 'AaBBCC',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')

    fig.text(0.87 , 0.19, 'AABBCC',  fontsize=11,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')

    plt.xticks(xticks_pos, freqs,  ha='right', rotation=45)
    plt.xticks([0,1,2, 3, 4, 5, 6], [0,1,2, 3, 4, 5, 6], rotation='horizontal')
    plt.xlabel('Number of height increasing alleles', fontsize = 18)
    plt.ylabel('Phenotype frequency', fontsize = 18)
    fig.tight_layout()
    fig_name = mydir + 'Figures/multiLocusTraits/L3.png'
    fig.savefig(fig_name, bbox_inches = None, pad_inches = 0.8, dpi = 600)
    plt.close()

    fig = plt.figure()
    #mu = 0
    #variance = 2
    #sigma = math.sqrt(variance)
    #x = np.linspace(-3, 3, 100)
    #plt.plot(x,mlab.normpdf(x, mu, sigma))
    np.random.seed(19680801)

    mu, sigma = 400, 50
    x = mu + sigma * np.random.randn(100000 )
    #x = 100*np.random.randn(100000)
    plt.xlabel('Number of height increasing alleles', fontsize = 18)
    plt.ylabel('Phenotype frequency', fontsize = 18)
    plt.hist(x, 100, fc='#87CEEB', histtype='bar', alpha=0.8, normed = 1)

    fig.tight_layout()
    fig_name = mydir + 'Figures/multiLocusTraits/L_inf.png'
    fig.savefig(fig_name, bbox_inches = None, pad_inches = 0.8, dpi = 600)
    plt.close()


def HWE():
    p = np.linspace(0, 1, num = 1000)
    p2 = p ** 2
    q = 1 - p
    q2 = q ** 2
    pq2 = 2 * p * q
    #fig = plt.figure()
    plt.subplot(111)

    plt.plot(p, p2, 'k-', lw = 4, c = '#87CEEB', label = r'$p^{2}$' )
    plt.plot(p, q2, 'k-', lw = 4, c = '#FFA500', label = r'$q^{2}$' )
    plt.plot(p, pq2, 'k-', lw = 4, c = '#FF6347', label = r'$2pq$' )
    plt.xlabel(r'$p$', fontsize = 18)
    plt.xlabel(r'$p$', fontsize = 18)
    plt.ylabel('Genotype frequency', fontsize = 18)
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
       ncol=3, mode="expand", borderaxespad=0.)

    fig_name = mydir + 'Figures/HWE.png'
    plt.savefig(fig_name, bbox_inches = None, pad_inches = 0.8, dpi = 600)
    plt.close()

def buriFig():
    IN = pd.read_csv(mydir + 'data/Buri/csv/Table13.csv', sep = ',', \
        header = 'infer', index_col = 0)
    to_keep = ['total_fixed_bw', '1', '2', '3', '4', '5', '6', '7', '8', \
            '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', \
            '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', \
            '31', 'total_fixed_bw75']
    IN_data = IN[to_keep]
    IN_data = IN_data[to_keep].div(IN_data.sum(axis=1), axis=0)
    fig, ax = plt.subplots()
    n = IN_data.shape[0]
    x = range(0, IN_data.shape[1])
    labels = ['0.0', '', '', '', '0.125', '', '', '', '0.25','', '', '', \
                '0.375', '', '', '', '0.5', '', '', '', '0.625', '', '', '', \
                '0.75', '', '', '', '0.875', '', '', '', '1.0']
    plt.xticks(x, labels, rotation='horizontal')

    barcollection = plt.bar(x,IN_data.ix[0,:].values, color = '#87CEEB')
    #time_text = ax.text(0.85, 0.9, '', transform=ax.transAxes)
    time_text = ax.text(.8, .9, '', fontsize=15)

    ax.set_title('''Buri's ''' + r'$D. \, melanogaster$' + ' drift experiment', fontsize = 22)
    ax.set_xlabel('Frequency of ' + r'$bw^{75}$' + ' allele', fontsize = 20)
    ax.set_ylabel('Fraction of populations', fontsize = 20)
    ax.set_xlim(0, 33)
    def animate(i):
        y = IN_data.ix[i,:].values
        time_text.set_text('Generation ' + str(i))
        for j, b in enumerate(barcollection):
            b.set_height(y[j])
        return time_text,
    anim = animation.FuncAnimation(fig, animate,repeat=False,blit=False,frames=n,
                                 interval=100)
    movie_name = mydir + 'Figures/Buri.gif'
    anim.save(movie_name, fps=1, writer='imagemagick')

def bumpusHist():
    IN = pd.read_csv(mydir + 'data/bumpus_full.csv', sep = ',', \
        header = 'infer')
    alive = IN.loc[IN['Survival'] == 'Alive'].Weight
    dead = IN.loc[IN['Survival'] == 'Dead'].Weight
    fig = plt.figure()
    plt.hist(alive, 15, fc='#87CEEB', histtype='bar', alpha=0.5, normed = False, label = 'Survived')
    plt.hist(dead, 15, fc='#FF6347', histtype='bar', alpha=0.5, normed = False, label = 'Died')
    plt.xlabel('Weight (grams)', fontsize = 18)
    plt.ylabel('Number of birds', fontsize = 18)
    plt.ylim(0, 20)
    plt.legend(loc = 'upper left')
    fig.tight_layout()
    fig_name = mydir + 'Figures/selection/bumpusHist.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


#driftFig(N = 100, n1 = 50,  s = 0, reps = 50)
#driftFig(N = 1000, n1 = 500,  s = 0, reps = 50)
#driftFig(N = 1000, n1 = 1,  s = 0.1, reps = 50, x_range = 150)
#lossHetero()
#multiLocusTraits(L = 2)
#galtonRegress()

bumpusHist()
