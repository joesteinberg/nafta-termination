########################################################################
# imports

import os
import re
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

colors=['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
hatches=[None,'x','/','.','+']
#colors=['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854']
import matplotlib.pyplot as plt

df = pd.read_pickle('output/wiod_trd.pik')

sectors=['A','R','T','M','S']
snames=['Agr.','Rsrcs.','Trans.','Mfg.','Srvcs.']

countries=['USA','CAN','MEX']
cnames=['USA','Canada','Mexico']

partners = countries + ['ROW']
pnames = cnames + ['ROW']

##########################################################################
# a little preprocessing
df=df[df.partner.isin(['USA','CAN','MEX','ROW'])]
df=df[df.sector!='TOT']
#df=df[df.year.isin(range(2010,2015))]
df=df[df.year==2014]

df['trd']=df.ex+df.im
df['trd_M']=df.ex_M+df.im_M
df['trd2'] = df.trd
df['va2'] = df.VA

cols = ['trd','trd2','trd_M','tb','va2']
for col in cols:
    if col in ['trd','trd_M','im_M','ex_M']:
        df[col]=100*df[col]/df.VA
    else:
        df[col]=100*df[col]/df.GDP

sum_s = df[df.partner!='ROW'].groupby(['region','sector','year'])['trd','trd2','trd_M','tb'].sum().reset_index()
avg_s = sum_s.groupby(['region','sector'])['trd','trd2','trd_M','tb'].mean().reset_index()
sum_s2 = df.groupby(['region','sector','year'])['trd','trd2','trd_M','tb'].sum().reset_index()
avg_s2 = sum_s2.groupby(['region','sector'])['trd','trd2','trd_M','tb'].mean().reset_index()

sum_p = df.groupby(['region','partner','year'])[cols].sum().reset_index()
avg_p = sum_p.groupby(['region','partner'])[cols].mean().reset_index()

avg_va = df.groupby(['region','sector'])['va2'].mean().reset_index()

##########################################################################
# figure 1: trade/GDP by partner

fig,axes=plt.subplots(1,1,figsize=(5,3.5))

data=[np.zeros(3) for p in partners]
inds=range(len(partners))

for i in inds:
    p=partners[i]
    for j in inds[0:3]:
        c=countries[j]
        tmp=avg_p.trd2[np.logical_and(avg_p.partner==p,avg_p.region==c)]
        if(len(tmp)>0):
            data[i][j]=tmp.values[0]

p1=axes.bar(inds[0:3],data[0],
            color=colors[0],edgecolor='black',linewidth=1,hatch=hatches[0],align='center',width=0.5)

p2=axes.bar(inds[0:3],data[1],bottom=data[0],
            color=colors[1],edgecolor='black',linewidth=1,hatch=4*hatches[1],align='center',width=0.5)

p3=axes.bar(inds[0:3],data[2],bottom=data[0]+data[1],
            color=colors[2],edgecolor='black',linewidth=1,hatch=4*hatches[2],align='center',width=0.5)

p4=axes.bar(inds[0:3],data[3],bottom=data[0]+data[1]+data[2],
            color=colors[3],edgecolor='black',linewidth=1,hatch=4*hatches[3],align='center',width=0.5)

axes.set_xticks(inds[0:4])
axes.set_xticklabels(countries)
axes.set_ylim(0,70)
axes.set_xlim(-0.5,2.5)
axes.set_ylabel('percent GDP')
#axes.set_title('(a) By partner',y=1.04,size=10)
axes.legend([p1,p2,p3,p4],pnames,loc='upper left',prop={'size':8},ncol=1)
plt.savefig('output/fig1_bilateral_trade.pdf',bbox='tight')

##########################################################################
# figure 2: sectoral trade relative to VA

width=0.25
fig,axes=plt.subplots(2,2,figsize=(7,7),sharex='row',sharey=False)

# (a) sectoral gross trade/GDP
inds=np.arange(len(sectors))

for i in range(len(countries)):
    data=np.zeros(len(inds))
    for j in range(len(sectors)):
        data[j]=avg_s.trd[np.logical_and(avg_s.region==countries[i],
                                         avg_s.sector==sectors[j])].values[0]          
    if i==0:
        p1=axes[0,0].bar(inds,data,align='edge',width=width,color=colors[0],hatch=hatches[0],linewidth=1)
    elif i==1:
        p2=axes[0,0].bar(inds+width,data,align='edge',width=width,color=colors[1],hatch=4*hatches[1],linewidth=1)
    elif i==2:
        p3=axes[0,0].bar(inds+2*width,data,align='edge',width=width,color=colors[2],hatch=4*hatches[2],linewidth=1)

axes[0,0].set_ylabel('percent sectoral value added')
axes[0,0].set_title('(a) Total NAFTA trade',y=1.04,size=10)
axes[0,0].set_xticks(inds+width*1.5)
axes[0,0].set_xticklabels(snames)
axes[0,0].set_xlim(-0.25,5.0)
axes[0,0].legend([p1,p2,p3],countries,loc='upper left',prop={'size':8})

# (c) sectoral intermediate trade/GDP
inds=np.arange(len(sectors))

for i in range(len(countries)):
    data=np.zeros(len(inds))
    for j in range(len(sectors)):
        data[j]=avg_s.trd_M[np.logical_and(avg_s.region==countries[i],
                                           avg_s.sector==sectors[j])].values[0]
          
    if i==0:
        p1=axes[0,1].bar(inds,data,align='edge',width=width,color=colors[0],hatch=hatches[0],linewidth=1)
    elif i==1:
        p2=axes[0,1].bar(inds+width,data,align='edge',width=width,color=colors[1],hatch=4*hatches[1],linewidth=1)
    elif i==2:
        p3=axes[0,1].bar(inds+2*width,data,align='edge',width=width,color=colors[2],hatch=4*hatches[2],linewidth=1)

axes[0,1].set_title('(b) Intermediate NAFTA trade',y=1.04,size=10)
#axes[0,1].set_xticks(inds+width*1.5)
#axes[0,1].set_xticklabels(snames)
#axes[0,1].set_xlim(-0.25,5.0)

#fig.subplots_adjust(hspace=0.0,wspace=0.15)
#plt.savefig('output/fig2_sectoral_trade_openness.pdf',bbox='tight')
#plt.clf()


width=0.25
#fig,axes=plt.subplots(1,2,figsize=(7,3.5),sharex=True,sharey=False)

# (a) sectoral VA/GDP
data=[np.zeros(3) for s in sectors]
inds=range(len(countries))

for i in range(len(sectors)):
    s=sectors[i]
    for j in inds:
        c=countries[j]
        tmp=avg_va.va2[np.logical_and(avg_va.sector==s,avg_va.region==c)]
        if(len(tmp)>0):
            data[i][j]=tmp.values[0]

p1=axes[1,0].bar(inds,data[0],
                 color=colors[0],hatch=hatches[0],linewidth=1,align='center',width=0.5)

p2=axes[1,0].bar(inds,data[1],bottom=data[0],
                 color=colors[1],hatch=4*hatches[1],linewidth=1,align='center',width=0.5)

p3=axes[1,0].bar(inds,data[2],bottom=data[0]+data[1],
                 color=colors[2],hatch=4*hatches[2],linewidth=1,align='center',width=0.5)

p4=axes[1,0].bar(inds,data[3],bottom=data[0]+data[1]+data[2],
                 color=colors[3],hatch=4*hatches[3],linewidth=1,align='center',width=0.5)

p5=axes[1,0].bar(inds,data[4],bottom=data[0]+data[1]+data[2]+data[3],
                 color=colors[4],hatch=4*hatches[4],linewidth=1,align='center',width=0.5)

axes[1,0].set_ylabel('percent GDP')
axes[1,0].set_title('(c) Value added shares',y=1.04,size=10)
axes[1,0].set_xticks(inds)
axes[1,0].set_xticklabels(countries)
axes[1,0].set_xlim(-0.5,2.5)
axes[1,0].set_ylim(0,140)
axes[1,0].legend([p1,p2,p3,p4,p5],snames,loc='upper left',prop={'size':8},ncol=2)

# (d) sectoral trade/GDP
data=[np.zeros(3) for s in sectors]
inds=range(len(countries))

for i in range(len(sectors)):
    s=sectors[i]
    for j in inds:
        c=countries[j]
        tmp=avg_s.trd2[np.logical_and(avg_s.sector==s,avg_s.region==c)]
        if(len(tmp)>0):
            data[i][j]=tmp.values[0]

p1=axes[1,1].bar(inds,data[0],
                 color=colors[0],hatch=hatches[0],linewidth=1,align='center',width=0.5)

p2=axes[1,1].bar(inds,data[1],bottom=data[0],
                 color=colors[1],hatch=4*hatches[1],linewidth=1,align='center',width=0.5)

p3=axes[1,1].bar(inds,data[2],bottom=data[0]+data[1],
                 color=colors[2],hatch=4*hatches[2],linewidth=1,align='center',width=0.5)

p4=axes[1,1].bar(inds,data[3],bottom=data[0]+data[1]+data[2],
                 color=colors[3],hatch=4*hatches[3],linewidth=1,align='center',width=0.5)

p5=axes[1,1].bar(inds,data[4],bottom=data[0]+data[1]+data[2]+data[3],
                 color=colors[4],hatch=4*hatches[4],linewidth=1,align='center',width=0.5)

axes[1,1].set_title('(d) Contributions to NAFTA trade',y=1.04,size=10)

fig.subplots_adjust(hspace=0.3,wspace=0.2)
plt.savefig('output/fig2_sectoral_trade.pdf',bbox='tight')
plt.clf()

# ##########################################################################
# # figure 3: net trade

vmax = lambda v: [max(x,0.0) for x in v]
vmin = lambda v: [min(x,0.0) for x in v]

def which_bottom(bottom_p, bottom_n, data):
    bottom=np.zeros(len(data))
    for i in range(len(data)):
        bottom[i]=bottom_p[i] if data[i]>0.0 else bottom_n[i]
    return bottom

fig,axes=plt.subplots(1,2,figsize=(7,3.5),sharex=True,sharey=True)

# (a) by partner
data=[np.zeros(3) for p in countries]
inds=range(len(countries))

for i in inds:
    p=partners[i]
    for j in inds[0:3]:
        c=countries[j]
        tmp=avg_p.tb[np.logical_and(avg_p.partner==p,avg_p.region==c)]
        if(len(tmp)>0):
            data[i][j]=tmp.values[0]

bottom_pos=np.zeros(3)
bottom_neg=np.zeros(3)
p1=axes[0].bar(inds[0:3],data[0],
               color=colors[0],hatch=hatches[0],linewidth=1,align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[0])
bottom_neg=bottom_neg+vmin(data[0])

bottom=which_bottom(bottom_pos,bottom_neg,data[1])
p2=axes[0].bar(inds[0:3],data[1],bottom=bottom,
               color=colors[1],hatch=4*hatches[1],linewidth=1,align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[1])
bottom_neg=bottom_neg+vmin(data[1])

bottom=which_bottom(bottom_pos,bottom_neg,data[2])
p3=axes[0].bar(inds[0:3],data[2],bottom=bottom,
               color=colors[2],hatch=4*hatches[2],linewidth=1,align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[2])
bottom_neg=bottom_neg+vmin(data[2])

axes[0].plot(range(-1,4),np.zeros(len(range(-1,4))),color='black',linestyle='-')
axes[0].set_xticks(inds)
axes[0].set_xticklabels(countries)
axes[0].set_yticks([-3,0,3,6,9,12])
axes[0].set_ylim(-3,12)
axes[0].set_xlim(-0.5,2.5)
axes[0].set_ylabel('percent GDP')
axes[0].set_title('(a) By partner',y=1.04,size=10)
axes[0].legend([p1,p2,p3],cnames,loc='upper left',prop={'size':8},ncol=1)

# (b) by sector
data=[np.zeros(3) for s in sectors]
inds=range(len(countries))

for i in range(len(sectors)):
    s=sectors[i]
    for j in inds:
        c=countries[j]
        tmp=avg_s.tb[np.logical_and(avg_s.sector==s,avg_s.region==c)]
        if(len(tmp)>0):
            data[i][j]=tmp.values[0]

bottom_pos=np.zeros(3)
bottom_neg=np.zeros(3)
p1=axes[1].bar(inds,data[0],color=colors[0],hatch=hatches[0],align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[0])
bottom_neg=bottom_neg+vmin(data[0])

bottom=which_bottom(bottom_pos,bottom_neg,data[1])
p2=axes[1].bar(inds,data[1],bottom=bottom,
               color=colors[1],hatch=4*hatches[1],linewidth=1,align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[1])
bottom_neg=bottom_neg+vmin(data[1])

bottom=which_bottom(bottom_pos,bottom_neg,data[2])
p3=axes[1].bar(inds,data[2],bottom=bottom,
               color=colors[2],hatch=4*hatches[2],linewidth=1,align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[2])
bottom_neg=bottom_neg+vmin(data[2])

bottom=which_bottom(bottom_pos,bottom_neg,data[3])
p4=axes[1].bar(inds,data[3],bottom=bottom,
               color=colors[3],hatch=4*hatches[3],linewidth=1,align='center',width=0.5)
bottom_pos=bottom_pos+vmax(data[3])
bottom_neg=bottom_neg+vmin(data[3])

bottom=which_bottom(bottom_pos,bottom_neg,data[4])
p5=axes[1].bar(inds,data[4],bottom=bottom,
               color=colors[4],hatch=4*hatches[4],linewidth=1,align='center',width=0.5)

axes[1].plot(range(-1,4),np.zeros(len(range(-1,4))),color='black',linestyle='-')
axes[1].set_title('(b) By sector',y=1.04,size=10)
axes[1].legend([p1,p2,p3,p4,p5],snames,loc='upper left',prop={'size':8},ncol=2)

fig.subplots_adjust(hspace=0.05,wspace=0.05)
plt.savefig('output/fig3_trade_balances.pdf',bbox='tight')
plt.clf()





plt.close('all')


