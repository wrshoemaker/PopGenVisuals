import popGenSims as pgs
import  matplotlib.pyplot as plt

def driftFig(N = 100, n1 = 50,  s= 0, reps = 5, mid_N = True):
    if mid_N == True:
        n1 = N / 2
    sims = pgs.multipleDriftSims(n1 = n1, N = N, s = s, reps = reps)
    for sim in sims:
        plt.plot(sim[0], sim[1], c = 'b')
    plt.show()
