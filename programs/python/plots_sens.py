########################################################################
# imports

import numpy as np
import pandas as pd

import matplotlib as mpl
mpl.use('Agg')
mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':16})
mpl.rc('font',size=16)
mpl.rc('text', usetex=True)
mpl.rc('lines',linewidth=1.5)
mpl.rc('savefig',bbox='tight')
mpl.rc('savefig',format='pdf')
mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':10})
mpl.rc('font',size=10)


hatches=[None,'x','/','.','+']
colors=['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
dashes = [(None,None),(10,2),(3,1)]
#colors=['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854']
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

sectors=['A','R','T','M','S']
snames=['Agr.','Rsrcs.','Trans.','Mfg.','Srvcs.']

countries=['USA','CAN','MEX']
cnames=['USA','Canada','Mexico']

partners = countries + ['ROW']
pnames = cnames + ['ROW']

def autolabel(ax,rects,sign):
    for rect in rects:
        height=sign*rect.get_height()
        ax.annotate('%0.2f' % float(height),
                    (rect.get_x()+rect.get_width()/2.0,height),
                    ha='center',
                    va='center',
                    xytext=(0,sign*20),
                    textcoords='offset points')


vmax = lambda v: [max(x,0.0) for x in v]
vmin = lambda v: [min(x,0.0) for x in v]

def which_bottom(bottom_p, bottom_n, data):
    bottom=np.zeros(len(data))
    for i in range(len(data)):
        bottom[i]=bottom_p[i] if data[i]>0.0 else bottom_n[i]
    return bottom

##########################################################################
# a little model results preprocessing

flag=2

suff=''
fignum=0

if flag==1:
    suff='_no_trd_adj_cost'
    fignum=7
elif flag==2:
    suff='_noio'
    fignum=8
else:
    raise

models_usa = []
models_usa.append(pd.read_csv('../c/output/vars0_usa'+suff+'.csv'))
models_usa.append(pd.read_csv('../c/output/vars1_usa'+suff+'.csv'))

models_can = []
models_can.append(pd.read_csv('../c/output/vars0_can'+suff+'.csv'))
models_can.append(pd.read_csv('../c/output/vars1_can'+suff+'.csv'))

models_mex = []
models_mex.append(pd.read_csv('../c/output/vars0_mex'+suff+'.csv'))
models_mex.append(pd.read_csv('../c/output/vars1_mex'+suff+'.csv'))

#models={'USA':models_usa,'CAN':models_can,'MEX':models_mex}
models=[models_usa,models_can,models_mex]
mperiods = models_usa[0].period

for m in models[0]:

    m['nx']=m['nx1']+m['nx2']
    m['trd_nafta']=m['ex1']+m['ex2']+m['im1']+m['im2']
    m['im_nafta']=m['im1']+m['im2']

    for s in [0,1,2,3,4]:
        m['exs'+str(s)+'_nafta']=m['exs'+str(s)+'-1']+m['exs'+str(s)+'-2']
        m['ims'+str(s)+'_nafta']=m['ims'+str(s)+'-1']+m['ims'+str(s)+'-2']
        m['nxs'+str(s)+'_nafta']=m['nxs'+str(s)+'-1']+m['nxs'+str(s)+'-2']
        m['imsm'+str(s)+'_nafta']=m['imsm'+str(s)+'-1']+m['imsm'+str(s)+'-2']
        #        m['exsm'+str(s)+'_nafta']=m['exsm'+str(s)+'-1']+m['exsm'+str(s)+'-2'
        m['nxs'+str(s)]=m['nxs'+str(s)+'_nafta']+m['nxs'+str(s)+'-3']
        m['exs'+str(s)]=m['exs'+str(s)+'_nafta']+m['exs'+str(s)+'-3']
        m['ims'+str(s)]=m['ims'+str(s)+'_nafta']+m['ims'+str(s)+'-3']

for m in models[1]:

    m['nx']=m['nx0']+m['nx2']
    m['trd_nafta']=m['ex0']+m['ex2']+m['im0']+m['im2']
    m['im_nafta']=m['im0']+m['im2']

    for s in [0,1,2,3,4]:
        m['exs'+str(s)+'_nafta']=m['exs'+str(s)+'-0']+m['exs'+str(s)+'-2']
        m['ims'+str(s)+'_nafta']=m['ims'+str(s)+'-0']+m['ims'+str(s)+'-2']
        m['nxs'+str(s)+'_nafta']=m['nxs'+str(s)+'-0']+m['nxs'+str(s)+'-2']
        m['imsm'+str(s)+'_nafta']=m['imsm'+str(s)+'-0']+m['imsm'+str(s)+'-2']
        #        m['exsm'+str(s)+'_nafta']=m['exsm'+str(s)+'-0']+m['exsm'+str(s)+'-2']
        m['nxs'+str(s)]=m['nxs'+str(s)+'_nafta']+m['nxs'+str(s)+'-3']
        m['exs'+str(s)]=m['exs'+str(s)+'_nafta']+m['exs'+str(s)+'-3']
        m['ims'+str(s)]=m['ims'+str(s)+'_nafta']+m['ims'+str(s)+'-3']
    
for m in models[2]:

    m['nx']=m['nx0']+m['nx1']
    m['trd_nafta']=m['ex0']+m['ex1']+m['im0']+m['im1']
    m['im_nafta']=m['im0']+m['im1']

    for s in [0,1,2,3,4]:
        m['exs'+str(s)+'_nafta']=m['exs'+str(s)+'-0']+m['exs'+str(s)+'-1']
        m['ims'+str(s)+'_nafta']=m['ims'+str(s)+'-0']+m['ims'+str(s)+'-1']
        m['nxs'+str(s)+'_nafta']=m['nxs'+str(s)+'-0']+m['nxs'+str(s)+'-1']
        m['imsm'+str(s)+'_nafta']=m['imsm'+str(s)+'-0']+m['imsm'+str(s)+'-1']
        #        m['exsm'+str(s)+'_nafta']=m['exsm'+str(s)+'-0']+m['exsm'+str(s)+'-1']
        m['nxs'+str(s)]=m['nxs'+str(s)+'_nafta']+m['nxs'+str(s)+'-3']
        m['exs'+str(s)]=m['exs'+str(s)+'_nafta']+m['exs'+str(s)+'-3']
        m['ims'+str(s)]=m['ims'+str(s)+'_nafta']+m['ims'+str(s)+'-3']


##########################################################################
# fig 6: dynamics

fig,axes=plt.subplots(3,2,figsize=(7,9),sharey=False,sharex=True)

for i in range(len(countries)):

    data = 100*(models[i][1]['im_nafta']/models[i][0]['im_nafta']-1.0)
    axes[0,0].plot(mperiods,data,color=colors[i],dashes=dashes[i])

    data = models[i][1]['te']
    axes[0,1].plot(mperiods,data,color=colors[i],dashes=dashes[i])
        
axes[0,0].set_xlim(0,40)
axes[0,0].set_title('(a) NAFTA imports',y=1.04,size=10)
axes[0,0].set_ylabel('pct. change')
axes[0,1].set_title('(b) Trade elasticity (log $\Delta$ imports/log $\Delta$ tariffs)',y=1.04,size=10)
axes[0,0].legend(cnames,loc='best',prop={'size':8},handlelength=3)

cols=['rgdp','c','ii','nx']
r=[1,1,2,2]
c=[0,1,0,1]
for i in range(len(cols)):
    col=cols[i]
    for j in range(len(countries)):
        if(col=='nx'):
            data=100*(models[j][1][col]/models[j][1].ngdp-models[j][0][col]/models[j][0].ngdp)
        else:
            data = 100*(models[j][1][col]/models[j][0][col]-1.0)
        axes[r[i],c[i]].plot(mperiods,data,color=colors[j],dashes=dashes[j])
        
    if():
        row=1

axes[2,0].set_xlabel('Period since termination')
axes[2,1].set_xlabel('Period since termination')
axes[1,0].set_xlim(0,40)
axes[1,0].set_title('(c) GDP',y=1.04,size=10)
axes[1,1].set_title('(d) Consumption',y=1.04,size=10)
axes[2,0].set_title('(e) Investment',y=1.04,size=10)
axes[2,1].set_title('(f) Net exports/GDP (p.p. chg.)',y=1.04,size=10)

axes[1,0].set_ylabel('pct. change')
axes[2,0].set_ylabel('pct. change')

axes[0,0].set_ylim(-16,0)
axes[0,1].set_ylim(0,12)
axes[1,0].set_ylim(-0.5,0.1)
axes[1,1].set_ylim(-0.3,0.0)
axes[2,0].set_ylim(-0.8,0.1)
axes[2,1].set_ylim(-0.2,0.5)


#axes[0,0].set_ylim(-16,0)
#axes[0,1].set_ylim(0,12)
#axes[1,0].set_ylim(-0.25,0.05)
#axes[1,1].set_ylim(-0.1,0.15)
#axes[2,0].set_ylim(-1,0.1)
#axes[2,1].set_ylim(-0.2,0.4)

fig.subplots_adjust(hspace=0.2,wspace=0.2)
plt.savefig('output/fig'+str(fignum)+'_trd_macro_dyn'+suff+'.pdf',bbox='tight')
plt.clf()
