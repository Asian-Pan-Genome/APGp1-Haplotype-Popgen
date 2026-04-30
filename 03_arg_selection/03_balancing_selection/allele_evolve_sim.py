from math import comb
import matplotlib.pyplot as plt
# import dadi
import numpy as np
from scipy.stats import binom
import os

import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def WF_transition_matrix(N):
    """
    Wright-Fisher transition matrix
    N+1 * N+1
    """
    P = np.zeros((N + 1, N + 1))
    
    for i in range(N + 1):
        p = i / N  # allele frequency
        P[i, :] = binom.pmf(np.arange(N + 1), N, p)
    
    return P

def WF_evolve(v, N, M):
    """
    Evolve one Wright-Fisher generation without explicitly forming
    the full transition matrix.
    """
    assert len(v) == N + 1
    
    w = np.zeros(M + 1)
    js = np.arange(M + 1)

    for i in range(N + 1):
        if v[i] == 0:
            continue
        
        p = i / N
        w += v[i] * binom.pmf(js, M, p)

    return w


def simulate():
    demo=[]
    #estimate from whole-genome ARG using 
    rates=np.array([266357.85246899,  38731.97174498,  17879.94395414,  12223.07070205,
            7882.51985768,   8732.08027795,  34127.13096024])
    # rates=[20000,20000,20000,20000,20000,20000]
    generations=np.array([500,1000,1250,1500,2000,np.inf])
    demo=(rates,generations)


    def prepare_forward_sim(demo,start_generation):
        rates,generations=demo
        idx = np.searchsorted(generations, start_generation, side="left")
        Ns=[]
        for j in range(idx+1):
            Ns.append(int(rates[j]))
        k=0
        all_sim=[]
        for j in range(idx+1):
            if j<idx:
                m=generations[j]-k
            else:
                m=start_generation-k
            all_sim.insert(0,[int(m),(Ns[j],Ns[j])])
            if j<idx:
                all_sim.insert(0,[1,(Ns[j+1],Ns[j])])
            k=generations[j]
        return Ns[::-1],all_sim

    p_thresholds=[0.1,0.13,0.15]
    all_curves=[[],[],[]]
    all_survivals=[[],[],[]]
    Ns,all_sim=prepare_forward_sim(demo,1500)
    #initial v Ns[0]*1, only v[1]=1
    v=np.zeros(Ns[0]+1)
    v[1]=1

    # conditional_tr=[]
    # survival=[]
    itr=0
    for l in all_sim:
        tn=l[1]
        if tn[0]!=tn[1]:
            #we use WF_evolve to get the next one without explicitly construction
            v=WF_evolve(v,tn[0],tn[1])
            for j in range(3):
                p_threshold=p_thresholds[j]
                idx_gt01=np.array([i for i in range(tn[1]+1)])/tn[1]>p_threshold
                all_curves[j].append(v[idx_gt01][:-1].sum()/v[1:-1].sum())
                all_survivals[j].append(v[1:-1].sum())
        else:
            P=WF_transition_matrix(tn[0])
            for i in range(int(l[0])):
                v=v@P
                itr+=1
                print(itr,"done")
                for j in range(3):
                    p_threshold=p_thresholds[j]
                    idx_gt01=np.array([i for i in range(tn[1]+1)])/tn[1]>p_threshold
                    all_curves[j].append(v[idx_gt01][:-1].sum()/v[1:-1].sum())
                    all_survivals[j].append(v[1:-1].sum())

    all_curves=np.array(all_curves)
    all_survivals=np.array(all_survivals)

    np.save(f"{SCRIPT_DIR}/discrete_sim_results/1500_all_curves.npy", np.array(all_curves))
    np.save(f"{SCRIPT_DIR}/discrete_sim_results/1500_all_survivals.npy", np.array(all_survivals))

def plot():
    p_thresholds=[0.1,0.13,0.15]
    all_curves=np.load(f"{SCRIPT_DIR}/discrete_sim_results/1500_all_curves.npy")

    plt.style.use("default")

    # plt.style.use("dark_background")

    fig, ax = plt.subplots(figsize=(9,5))
    colors = plt.cm.Blues(np.linspace(0.4, 0.9, len(p_thresholds)))

    for i in range(0,3):
        ax.plot(
            all_curves[i],
            label=f"f = {p_thresholds[i]}",
            linewidth=2,
            color=colors[i]
        )
    ymin, ymax = plt.ylim()
    yticks = np.arange(
        0,
        (np.ceil(ymax / 0.005)+1) * 0.005,
        0.005
    )
    ax.set_yticks(yticks)

    # Show horizontal grid only
    ax.grid(axis='y', which='major', linestyle='--', alpha=0.5)


    for x in [250, 500, 1000,1500]:
        plt.axvline(x=x, linestyle='--', linewidth=1, alpha=0.5, zorder=0,c="grey")
    ax.set_xticks([0,250, 500, 1000,1500])


    TICK_SIZE = 12
    LABEL_SIZE = 14
    LEGEND_SIZE = 12
    TITLE_SIZE = 15

    ax.tick_params(axis='both', labelsize=TICK_SIZE)
    ax.set_xlabel("Generation", fontsize=LABEL_SIZE)
    ax.set_ylabel("P(p(T)>f | segregating)", fontsize=LABEL_SIZE)
    # ax.set_title("Probability of Reaching Frequency f Over Generations",
    #              fontsize=TITLE_SIZE, pad=10)
    ax.legend(fontsize=LEGEND_SIZE, frameon=False,loc="upper left")

    plt.savefig(
        f"{SCRIPT_DIR}/probability_curve.pdf",
        bbox_inches="tight"
    )

    plt.show()

def main():
    simulate()
    plot()

if __name__ == "__main__":
    main()