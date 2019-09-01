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

suff=''
path = 'c'

models_usa = []
models_usa.append(pd.read_csv('../'+path+'/output/vars0_usa'+suff+'.csv'))
models_usa.append(pd.read_csv('../'+path+'/output/vars1_usa'+suff+'.csv'))

models_can = []
models_can.append(pd.read_csv('../'+path+'/output/vars0_can'+suff+'.csv'))
models_can.append(pd.read_csv('../'+path+'/output/vars1_can'+suff+'.csv'))

models_mex = []
models_mex.append(pd.read_csv('../'+path+'/output/vars0_mex'+suff+'.csv'))
models_mex.append(pd.read_csv('../'+path+'/output/vars1_mex'+suff+'.csv'))

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
# figure 5: LR trade

fig,axes=plt.subplots(3,2,figsize=(7,9),sharex=True,sharey=False)

# (a) gross trade by partner
data=[np.zeros(3) for p in partners]
inds=range(len(partners))

for i in inds:
    p=partners[i]
    for j in inds[0:3]:
        if(j!=i):
            c=countries[j]
            tmp0=(models[j][0]['im'+str(i)][50]+models[j][0]['ex'+str(i)][50])/models[j][0].ngdp[50]
            tmp1=(models[j][1]['im'+str(i)][50]+models[j][1]['ex'+str(i)][50])/models[j][1].ngdp[50]
            data[i][j]=100.0*(tmp1-tmp0)

axes[0,0].axhline(0.0,color='black',linewidth=1,linestyle='-')

bottom_pos=np.zeros(3)
bottom_neg=np.zeros(3)
p1=axes[0,0].bar(inds[0:3],data[0],
                 color=colors[0],edgecolor='black',linewidth=1,hatch=hatches[0],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[0])
bottom_neg=bottom_neg+vmin(data[0])

bottom=which_bottom(bottom_pos,bottom_neg,data[1])
p2=axes[0,0].bar(inds[0:3],data[1],bottom=bottom,
                 color=colors[1],edgecolor='black',linewidth=1,hatch=4*hatches[1],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[1])
bottom_neg=bottom_neg+vmin(data[1])

bottom=which_bottom(bottom_pos,bottom_neg,data[2])
p3=axes[0,0].bar(inds[0:3],data[2],bottom=bottom,
                 color=colors[2],edgecolor='black',linewidth=1,hatch=4*hatches[2],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[2])
bottom_neg=bottom_neg+vmin(data[2])

bottom=which_bottom(bottom_pos,bottom_neg,data[3])
p4=axes[0,0].bar(inds[0:3],data[3],bottom=bottom,
                 color=colors[3],edgecolor='black',linewidth=1,hatch=4*hatches[3],align='center',width=0.5)

axes[0,0].set_xticks(inds[0:4])
axes[0,0].set_xticklabels(countries)
#axes.set_ylim(0,70)
axes[0,0].set_xlim(-0.5,2.5)
axes[0,0].set_ylabel('pct. GDP, change')
axes[0,0].set_title('(a) Bilateral gross trade',y=1.04,size=10)
axes[0,0].legend([p1,p2,p3,p4],pnames,loc='lower left',prop={'size':8},ncol=1)

# (b) net trade by partner
data=[np.zeros(3) for p in partners]
inds=range(len(partners))

for i in inds:
    p=partners[i]
    for j in inds[0:3]:
        if(j!=i):
            c=countries[j]
            tmp0=(models[j][0]['nx'+str(i)][50])/models[j][0].ngdp[50]
            tmp1=(models[j][1]['nx'+str(i)][50])/models[j][1].ngdp[50]
            data[i][j]=100.0*(tmp1-tmp0)

axes[0,1].axhline(0.0,color='black',linewidth=1,linestyle='-')

bottom_pos=np.zeros(3)
bottom_neg=np.zeros(3)
p1=axes[0,1].bar(inds[0:3],data[0],
                 color=colors[0],edgecolor='black',linewidth=1,hatch=hatches[0],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[0])
bottom_neg=bottom_neg+vmin(data[0])

bottom=which_bottom(bottom_pos,bottom_neg,data[1])
p2=axes[0,1].bar(inds[0:3],data[1],bottom=bottom,
                 color=colors[1],edgecolor='black',linewidth=1,hatch=4*hatches[1],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[1])
bottom_neg=bottom_neg+vmin(data[1])

bottom=which_bottom(bottom_pos,bottom_neg,data[2])
p3=axes[0,1].bar(inds[0:3],data[2],bottom=bottom,
                 color=colors[2],edgecolor='black',linewidth=1,hatch=4*hatches[2],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[2])
bottom_neg=bottom_neg+vmin(data[2])

bottom=which_bottom(bottom_pos,bottom_neg,data[3])
p4=axes[0,1].bar(inds[0:3],data[3],bottom=bottom,
                 color=colors[3],edgecolor='black',linewidth=1,hatch=4*hatches[3],align='center',width=0.5)

axes[0,1].set_title('(b) Bilateral trade imbalances',y=1.04,size=10)
axes[0,1].set_ylim(-0.5,0.5)
#axes[0,1].legend([p1,p2,p3,p4],pnames,loc='lower left',prop={'size':8},ncol=1)


# (c) gross NAFTA trade by sector
data=[np.zeros(3) for s in sectors]
inds=range(len(sectors))

for i in inds:
    for j in inds[0:3]:
        tmp0=(models[j][0]['ims'+str(i)+'_nafta'][50]+models[j][0]['exs'+str(i)+'_nafta'][50])/models[j][0].ngdp[50]
        tmp1=(models[j][1]['ims'+str(i)+'_nafta'][50]+models[j][1]['exs'+str(i)+'_nafta'][50])/models[j][1].ngdp[50]
        data[i][j]=100.0*(tmp1-tmp0)

axes[1,0].axhline(0.0,color='black',linewidth=1,linestyle='-')

bottom_pos=np.zeros(3)
bottom_neg=np.zeros(3)
p1=axes[1,0].bar(inds[0:3],data[0],
                 color=colors[0],edgecolor='black',linewidth=1,hatch=hatches[0],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[0])
bottom_neg=bottom_neg+vmin(data[0])

bottom=which_bottom(bottom_pos,bottom_neg,data[1])
p2=axes[1,0].bar(inds[0:3],data[1],bottom=bottom,
                 color=colors[1],edgecolor='black',linewidth=1,hatch=4*hatches[1],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[1])
bottom_neg=bottom_neg+vmin(data[1])

bottom=which_bottom(bottom_pos,bottom_neg,data[2])
p3=axes[1,0].bar(inds[0:3],data[2],bottom=bottom,
                 color=colors[2],edgecolor='black',linewidth=1,hatch=4*hatches[2],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[2])
bottom_neg=bottom_neg+vmin(data[2])

bottom=which_bottom(bottom_pos,bottom_neg,data[3])
p4=axes[1,0].bar(inds[0:3],data[3],bottom=bottom,
                 color=colors[3],edgecolor='black',linewidth=1,hatch=4*hatches[3],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[3])
bottom_neg=bottom_neg+vmin(data[3])

bottom=which_bottom(bottom_pos,bottom_neg,data[4])
p5=axes[1,0].bar(inds[0:3],data[4],bottom=bottom,
                 color=colors[4],edgecolor='black',linewidth=1,hatch=4*hatches[4],align='center',width=0.5)

axes[1,0].set_ylabel('pct. GDP, change')
axes[1,0].set_title('(c) NAFTA sectoral gross trade',y=1.04,size=10)
axes[1,0].legend([p1,p2,p3,p4,p5],snames,loc='lower left',prop={'size':8},ncol=1)

# (d) net ROW trade by sector
data=[np.zeros(3) for s in sectors]
inds=range(len(sectors))

for i in inds:
    for j in inds[0:3]:
        tmp0=(models[j][0]['nxs'+str(i)+'_nafta'][50])/models[j][0].ngdp[50]
        tmp1=(models[j][1]['nxs'+str(i)+'_nafta'][50])/models[j][1].ngdp[50]
        data[i][j]=100.0*(tmp1-tmp0)

axes[1,1].axhline(0.0,color='black',linewidth=1,linestyle='-')

bottom_pos=np.zeros(3)
bottom_neg=np.zeros(3)
p1=axes[1,1].bar(inds[0:3],data[0],
                 color=colors[0],edgecolor='black',linewidth=1,hatch=hatches[0],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[0])
bottom_neg=bottom_neg+vmin(data[0])

bottom=which_bottom(bottom_pos,bottom_neg,data[1])
p2=axes[1,1].bar(inds[0:3],data[1],bottom=bottom,
                 color=colors[1],edgecolor='black',linewidth=1,hatch=4*hatches[1],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[1])
bottom_neg=bottom_neg+vmin(data[1])

bottom=which_bottom(bottom_pos,bottom_neg,data[2])
p3=axes[1,1].bar(inds[0:3],data[2],bottom=bottom,
                 color=colors[2],edgecolor='black',linewidth=1,hatch=4*hatches[2],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[2])
bottom_neg=bottom_neg+vmin(data[2])

bottom=which_bottom(bottom_pos,bottom_neg,data[3])
p4=axes[1,1].bar(inds[0:3],data[3],bottom=bottom,
                 color=colors[3],edgecolor='black',linewidth=1,hatch=4*hatches[3],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[3])
bottom_neg=bottom_neg+vmin(data[3])

bottom=which_bottom(bottom_pos,bottom_neg,data[4])
p5=axes[1,1].bar(inds[0:3],data[4],bottom=bottom,
                 color=colors[4],edgecolor='black',linewidth=1,hatch=4*hatches[4],align='center',width=0.5)

#axes[1,1].set_ylim(-0.4,0.6)
axes[1,1].set_title('(d) NAFTA sectoral trade imbalances',y=1.04,size=10)
#axes[1,1].legend([p1,p2,p3,p4,p5],snames,loc='lower left',prop={'size':8},ncol=1)

# (e) gross ROW trade by sector
data=[np.zeros(3) for s in sectors]
inds=range(len(sectors))

for i in inds:
    for j in inds[0:3]:
        tmp0=(models[j][0]['ims'+str(i)+'-3'][50]+
              models[j][0]['exs'+str(i)+'-3'][50])/models[j][0].ngdp[50]
        tmp1=(models[j][1]['ims'+str(i)+'-3'][50]+
              models[j][1]['exs'+str(i)+'-3'][50])/models[j][1].ngdp[50]
        data[i][j]=100.0*(tmp1-tmp0)

axes[2,0].axhline(0.0,color='black',linewidth=1,linestyle='-')

bottom_pos=np.zeros(3)
bottom_neg=np.zeros(3)
p1=axes[2,0].bar(inds[0:3],data[0],
                 color=colors[0],edgecolor='black',linewidth=1,hatch=hatches[0],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[0])
bottom_neg=bottom_neg+vmin(data[0])

bottom=which_bottom(bottom_pos,bottom_neg,data[1])
p2=axes[2,0].bar(inds[0:3],data[1],bottom=bottom,
                 color=colors[1],edgecolor='black',linewidth=1,hatch=4*hatches[1],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[1])
bottom_neg=bottom_neg+vmin(data[1])

bottom=which_bottom(bottom_pos,bottom_neg,data[2])
p3=axes[2,0].bar(inds[0:3],data[2],bottom=bottom,
                 color=colors[2],edgecolor='black',linewidth=1,hatch=4*hatches[2],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[2])
bottom_neg=bottom_neg+vmin(data[2])

bottom=which_bottom(bottom_pos,bottom_neg,data[3])
p4=axes[2,0].bar(inds[0:3],data[3],bottom=bottom,
                 color=colors[3],edgecolor='black',linewidth=1,hatch=4*hatches[3],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[3])
bottom_neg=bottom_neg+vmin(data[3])

bottom=which_bottom(bottom_pos,bottom_neg,data[4])
p5=axes[2,0].bar(inds[0:3],data[4],bottom=bottom,
                 color=colors[4],edgecolor='black',linewidth=1,hatch=4*hatches[4],align='center',width=0.5)

axes[2,0].set_ylabel('pct. GDP, change')
#axes[2,0].set_ylim(-0.1,0.6)
axes[2,0].set_title('(e) ROW sectoral gross trade',y=1.04,size=10)
axes[2,0].legend([p1,p2,p3,p4,p5],snames,loc='upper left',prop={'size':8},ncol=1)

# (f) net ROW trade by sector
data=[np.zeros(3) for s in sectors]
inds=range(len(sectors))

for i in inds:
    for j in inds[0:3]:
        tmp0=(models[j][0]['nxs'+str(i)+'-3'][50])/models[j][0].ngdp[50]
        tmp1=(models[j][1]['nxs'+str(i)+'-3'][50])/models[j][1].ngdp[50]
        data[i][j]=100.0*(tmp1-tmp0)

axes[2,1].axhline(0.0,color='black',linewidth=1,linestyle='-')

bottom_pos=np.zeros(3)
bottom_neg=np.zeros(3)
p1=axes[2,1].bar(inds[0:3],data[0],
                 color=colors[0],edgecolor='black',linewidth=1,hatch=hatches[0],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[0])
bottom_neg=bottom_neg+vmin(data[0])

bottom=which_bottom(bottom_pos,bottom_neg,data[1])
p2=axes[2,1].bar(inds[0:3],data[1],bottom=bottom,
                 color=colors[1],edgecolor='black',linewidth=1,hatch=4*hatches[1],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[1])
bottom_neg=bottom_neg+vmin(data[1])

bottom=which_bottom(bottom_pos,bottom_neg,data[2])
p3=axes[2,1].bar(inds[0:3],data[2],bottom=bottom,
                 color=colors[2],edgecolor='black',linewidth=1,hatch=4*hatches[2],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[2])
bottom_neg=bottom_neg+vmin(data[2])

bottom=which_bottom(bottom_pos,bottom_neg,data[3])
p4=axes[2,1].bar(inds[0:3],data[3],bottom=bottom,
                 color=colors[3],edgecolor='black',linewidth=1,hatch=4*hatches[3],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[3])
bottom_neg=bottom_neg+vmin(data[3])

bottom=which_bottom(bottom_pos,bottom_neg,data[4])
p5=axes[2,1].bar(inds[0:3],data[4],bottom=bottom,
                 color=colors[4],edgecolor='black',linewidth=1,hatch=4*hatches[4],align='center',width=0.5)

#axes[2,1].set_ylim(-0.6,0.4)
axes[2,1].set_title('(f) ROW sectoral trade imbalances',y=1.04,size=10)
#axes[2,1].legend([p1,p2,p3,p4,p5],snames,loc='lower left',prop={'size':8},ncol=1)


fig.subplots_adjust(hspace=0.3,wspace=0.2)
plt.savefig('output/fig4_LR_trade.pdf',bbox='tight')

##########################################################################
# figure 6: LR sectoral reallocation

fig,axes=plt.subplots(1,2,figsize=(7,3.5),sharex=True,sharey=False)

# (e) gross ROW trade by sector
data=[np.zeros(3) for s in sectors]
inds=range(len(sectors))

for i in inds:
    for j in inds[0:3]:
        tmp0=(models[j][0]['va'+str(i)][50])/models[j][0].ngdp[50]
        tmp1=(models[j][1]['va'+str(i)][50])/models[j][0].ngdp[50]
        data[i][j]=100.0*(tmp1-tmp0)

axes[0].axhline(0.0,color='black',linewidth=1,linestyle='-')

bottom_pos=np.zeros(3)
bottom_neg=np.zeros(3)
p1=axes[0].bar(inds[0:3],data[0],
            color=colors[0],edgecolor='black',linewidth=1,hatch=hatches[0],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[0])
bottom_neg=bottom_neg+vmin(data[0])

bottom=which_bottom(bottom_pos,bottom_neg,data[1])
p2=axes[0].bar(inds[0:3],data[1],bottom=bottom,
            color=colors[1],edgecolor='black',linewidth=1,hatch=4*hatches[1],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[1])
bottom_neg=bottom_neg+vmin(data[1])

bottom=which_bottom(bottom_pos,bottom_neg,data[2])
p3=axes[0].bar(inds[0:3],data[2],bottom=bottom,
            color=colors[2],edgecolor='black',linewidth=1,hatch=4*hatches[2],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[2])
bottom_neg=bottom_neg+vmin(data[2])

bottom=which_bottom(bottom_pos,bottom_neg,data[3])
p4=axes[0].bar(inds[0:3],data[3],bottom=bottom,
            color=colors[3],edgecolor='black',linewidth=1,hatch=4*hatches[3],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[3])
bottom_neg=bottom_neg+vmin(data[3])

bottom=which_bottom(bottom_pos,bottom_neg,data[4])
p5=axes[0].bar(inds[0:3],data[4],bottom=bottom,
               color=colors[4],edgecolor='black',linewidth=1,hatch=4*hatches[4],align='center',width=0.5)

axes[0].set_xticks(inds[0:4])
axes[0].set_xticklabels(countries)
axes[0].set_xlim(-0.5,2.5)

axes[0].set_ylabel('pct. initial aggregate, change')
axes[0].set_ylim(-1,0.5)
axes[0].set_title('(a) Sectoral value added',y=1.04,size=10)
axes[0].legend([p1,p2,p3,p4,p5],snames,loc='lower left',prop={'size':8},ncol=1)

# (b) net ROW trade by sector
data=[np.zeros(3) for s in sectors]
inds=range(len(sectors))

for i in inds:
    for j in inds[0:3]:
        tmp0=(models[j][0]['c'+str(i)][50])/models[j][0].c[50]
        tmp1=(models[j][1]['c'+str(i)][50])/models[j][0].c[50]
        data[i][j]=100.0*(tmp1-tmp0)

axes[1].axhline(0.0,color='black',linewidth=1,linestyle='-')

bottom_pos=np.zeros(3)
bottom_neg=np.zeros(3)
p1=axes[1].bar(inds[0:3],data[0],
            color=colors[0],edgecolor='black',linewidth=1,hatch=hatches[0],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[0])
bottom_neg=bottom_neg+vmin(data[0])

bottom=which_bottom(bottom_pos,bottom_neg,data[1])
p2=axes[1].bar(inds[0:3],data[1],bottom=bottom,
            color=colors[1],edgecolor='black',linewidth=1,hatch=4*hatches[1],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[1])
bottom_neg=bottom_neg+vmin(data[1])

bottom=which_bottom(bottom_pos,bottom_neg,data[2])
p3=axes[1].bar(inds[0:3],data[2],bottom=bottom,
            color=colors[2],edgecolor='black',linewidth=1,hatch=4*hatches[2],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[2])
bottom_neg=bottom_neg+vmin(data[2])

bottom=which_bottom(bottom_pos,bottom_neg,data[3])
p4=axes[1].bar(inds[0:3],data[3],bottom=bottom,
            color=colors[3],edgecolor='black',linewidth=1,hatch=4*hatches[3],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[3])
bottom_neg=bottom_neg+vmin(data[3])

bottom=which_bottom(bottom_pos,bottom_neg,data[4])
p5=axes[1].bar(inds[0:3],data[4],bottom=bottom,
               color=colors[4],edgecolor='black',linewidth=1,hatch=4*hatches[4],align='center',width=0.5)

#axes[1].set_ylim(-0.6,0.4)
axes[1].set_title('(b) Sectoral consumption',y=1.04,size=10)
#axes[2,1].legend([p1,p2,p3,p4,p5],snames,loc='lower left',prop={'size':8},ncol=1)


fig.subplots_adjust(hspace=0.3,wspace=0.2)
plt.savefig('output/fig5_LR_realloc.pdf',bbox='tight')


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
#axes[0,1].set_ylabel('log(chg. imports)/log(chg. tariffs)')
axes[0,0].legend(cnames,loc='best',prop={'size':8},handlelength=3)

#fig.subplots_adjust(hspace=0.0,wspace=0.25)
#plt.savefig('output/fig6_trd_dyn.pdf',bbox='tight')
#plt.clf()

#fig,axes=plt.subplots(2,2,figsize=(7,6.5),sharey=False,sharex=True)

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
    #axes[r[i],c[i]].axhline(0.0,color='black',linewidth=1,alpha=0.5,linestyle=':')
        
    if():
        row=1

axes[2,0].set_xlabel('Period since termination')
axes[2,1].set_xlabel('Period since termination')
#axes[1,0].set_xlim(0,40)
axes[1,0].set_title('(c) GDP',y=1.04,size=10)
#axes[1,1].set_ylim(-0.1,0.15)
axes[1,1].set_title('(d) Consumption',y=1.04,size=10)
axes[2,0].set_title('(e) Investment',y=1.04,size=10)
axes[2,1].set_title('(f) Net exports/GDP (p.p. chg.)',y=1.04,size=10)

axes[1,0].set_ylabel('pct. change')
axes[2,0].set_ylabel('pct. change')

#axes[0,0].set_ylim(-16,0)
#axes[0,1].set_ylim(0,12)
#axes[1,0].set_ylim(-0.25,0.05)
#axes[1,1].set_ylim(-0.1,0.15)
#axes[2,0].set_ylim(-1,0.1)
#axes[2,1].set_ylim(-0.2,0.4)

#axes[0,0].legend(cnames,loc='best',prop={'size':8})

fig.subplots_adjust(hspace=0.2,wspace=0.2)
plt.savefig('output/fig6_trd_macro_dyn.pdf',bbox='tight')
plt.clf()


plt.clf()
plt.close('all')
