"""

Plot Figures 4C & 4D from:
Berecki, Bryson et al. SCN1A Gain of Function in Early Infantile Encephalopathy, Ann Neurol 2019

Script takes 'trace' or 'if' as argument to plot time-voltage trace or I-F relationship. Ie

python3 run.py --trace

"""

import sys,os,argparse,neuron
from neuron import h
import numpy as np
import matplotlib.pyplot as plt


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--plot', required=True)
    args,unknown = parser.parse_known_args()
    PLOT = str(args.plot)

    # load PV template
    h.load_file('PV_template.hoc')

    # plot time-voltage traces
    if PLOT=='trace':
        params={'Het':[0.5,0.3],'Hom':[1.0,0.175]} # specify proportion of 'mutated' channels and current input
        Traces(params)

    # plot IF curve
    if PLOT == 'if':
        inputs = np.linspace(0,0.15,40) # input current (nA)
        IF(inputs)

    print('done')

# GET TIME-VOLTAGE TRACES
def Traces(params):
    for P,p in params.items():

        # create cells
        PvM,PvWT=h.pv('morphologies','C210401C.asc'),h.pv('morphologies','C210401C.asc')

        # get T226m trace
        PvM = mut(PvM,p[0]) # 'mutate' Nav1.1
        tM,vM,apM= getTrace(PvM,p[1],180)
        fig,ax=plt.subplots(1,1,figsize=(5,2.0))
        ax.plot(tM,vM,'red')

        # get WT trace
        tWT,vWT,apWT = getTrace(PvWT,p[1],180)
        ax.plot(tWT,vWT,'black'), ax.legend(['T226M', 'WT'],fontsize='x-small',loc='upper right')
        ax.set_title(P,weight='demi'),ax.set_ylabel('mV',fontsize=8),ax.set_xlabel('ms')

        fig.tight_layout()
        fig.savefig('save/'+P+'.pdf')


# GET IF CURVE
def IF(inputs):

    # create cells
    [PvHet,PvHom,PvWT]=[h.pv('morphologies','C210401C.asc') for X in range(3)]

    # 'mutate' proportion of Nav11 channels
    PvHom = mut(PvHom,1.0)
    PvHet = mut(PvHet,0.5)

    # get WT IF
    apWT = getIF(inputs,PvWT)

    # get Het IF
    apHet = getIF(inputs,PvHet)

    # get Hom IF
    apHom = getIF(inputs,PvHom)

    # plot
    fig,ax=plt.subplots(1,2,figsize=(10,3.0))
    ax[0].plot(inputs*1000,apWT,'o-',color='black'),ax[0].plot(inputs*1000,apHet,'o-',color='red')
    ax[1].plot(inputs*1000,apWT,'o-',color='black'),ax[1].plot(inputs*1000,apHom,'o-',color='red')
    for AX in ax:
        AX.set_xlabel('pA',fontsize=9),AX.set_ylabel('Hz',fontsize=9)
        AX.legend(['WT','T226M'])
    ax[0].set_title('Het', weight='demi'),ax[1].set_title('Hom', weight='demi')
    fig.tight_layout()

    # save
    fig.savefig('save/if.pdf')


#  *** OTHER FUNCTIONS ***

# change Nav1.1 conductance within 'mutated' Nav11m channels
def mut(Pv,MUT):
    for sec in Pv.all:
        if 'Nav11' in [mech.name() for seg in sec.allseg() for mech in seg]:
            gNav = sec(0.5).gNav11bar_Nav11
        for seg in sec:
            for mech in seg:
                if mech.name()=='Nav11m':
                    seg.gNav11bar_Nav11m = MUT*gNav
                    seg.mh_Nav11m, seg.hh_Nav11m = -26.6, -60.2
                    seg.tmh_Nav11m,seg.thh_Nav11m = -40.0, -65.0
                if mech.name()=='Nav11':
                    seg.gNav11bar_Nav11 = (1.0 - MUT)*gNav
    return Pv

def Icl(sec,AMP,DUR):
    Ic = h.IClamp(sec(0.5))
    Ic.dur, Ic.delay, Ic.amp = DUR,20,AMP
    return Ic

def Vec(sec, VAR):
    v = h.Vector()
    if VAR=='V': v.record(sec(0.5)._ref_v)
    if VAR=='T':v.record(h._ref_t)
    if VAR=='APC':v = h.APCount(sec(0.5))
    return v

def getTrace(Pv,AMP,DUR):
    stim=Icl(Pv.soma[0],AMP,DUR)
    t,v,APc=Vec(Pv.soma[0],'T'),Vec(Pv.soma[0],'V'),Vec(Pv.soma[0],'APC')
    hRun(DUR+20)
    return t,v,APc

def getIF(inputs,Pv):
    aps = []
    for AMP in inputs: aps.append(getTrace(Pv,AMP,500)[2].n*2.0) # run over 500ms
    return aps

def hRun(T):
    h.tstop,h.celsius=T,24.0
    h.cvode_active(1)
    h.run()

if __name__=='__main__':
    main()
