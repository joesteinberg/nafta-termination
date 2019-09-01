#####################################################################################
#imports, small functions, etc.

import numpy as np
from scipy import signal
import pandas as pd

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':26})
mpl.rc('font',size=26)
mpl.rc('text', usetex=True)
mpl.rc('lines',linewidth=1.5)
mpl.rc('savefig',bbox='tight')
mpl.rc('savefig',format='pdf')

cnt = 0

def norm(x):
    return 100*x/x.iloc[0]

#################################################################################
# load the model results

models_usa = []
models_usa.append(pd.read_csv('../c/output/vars0_usa.csv'))
models_usa.append(pd.read_csv('../c/output/vars1_usa.csv'))

models_can = []
models_can.append(pd.read_csv('../c/output/vars0_can.csv'))
models_can.append(pd.read_csv('../c/output/vars1_can.csv'))

models_mex = []
models_mex.append(pd.read_csv('../c/output/vars0_mex.csv'))
models_mex.append(pd.read_csv('../c/output/vars1_mex.csv'))

models={'USA':models_usa,'CAN':models_can,'MEX':models_mex}
mperiods = models_usa[0].period
mseries=[]

for m in models['USA']:
    m['nx']=m.nx1+m.nx2+m.nx3

    m['ex']=m.rex1+m.rex2+m.rex3
    m['im']=m.rim1+m.rim2+m.rim3
    m['trd']=m.ex+m.im
    m['reer']=(m.rer1*(m.rex1+m.rim1) + m.rer2*(m.rex2+m.rim2) + m.rer3*(m.rex3+m.rim3))/(m.rex1+m.rim1+m.rex2+m.rim2+m.rex3+m.rim3)

    for s in [0,1,2,3,4]:
        m['exs'+str(s)]=m['rexs'+str(s)+'-1']+m['rexs'+str(s)+'-2']+m['rexs'+str(s)+'-3']
        m['ims'+str(s)]=m['rims'+str(s)+'-1']+m['rims'+str(s)+'-2']+m['rims'+str(s)+'-3']
        m['trds'+str(s)]=m['exs'+str(s)]+m['ims'+str(s)]
        m['nxs'+str(s)]=m['nxs'+str(s)+'-1']+m['nxs'+str(s)+'-2']+m['nxs'+str(s)+'-3']

    m['ex_nafta']=m.rex1+m.rex2
    m['im_nafta']=m.rim1+m.rim2
    m['trd_nafta']=m.ex_nafta+m.im_nafta
    m['nx_nafta']=m.nx1+m.nx2
    m['reer_nafta']=(m.rer1*(m.rex1+m.rim1) + m.rer2*(m.rex2+m.rim2))/(m.rex1+m.rim1+m.rex2+m.rim2)


    for s in [0,1,2,3,4]:
        m['exs'+str(s)+'_nafta']=m['rexs'+str(s)+'-1']+m['rexs'+str(s)+'-2']
        m['ims'+str(s)+'_nafta']=m['rims'+str(s)+'-1']+m['rims'+str(s)+'-2']
        m['trds'+str(s)+'_nafta']=m['exs'+str(s)+'_nafta']+m['ims'+str(s)+'_nafta']
        m['nxs'+str(s)+'_nafta']=m['nxs'+str(s)+'-1']+m['nxs'+str(s)+'-2']


    m['trd_row']=m.rex3+m.rim3

for m in models['CAN']:
    m['nx']=m.nx0+m.nx2+m.nx3

    m['ex']=m.rex0+m.rex2+m.rex3
    m['im']=m.rim0+m.rim2+m.rim3
    m['trd']=m.ex+m.im
    m['reer']=(m.rer0*(m.rex0+m.rim0) + m.rer2*(m.rex2+m.rim2) + m.rer3*(m.rex3+m.rim3))/(m.rex0+m.rim0+m.rex2+m.rim2+m.rex3+m.rim3)

    for s in [0,1,2,3,4]:
        m['exs'+str(s)]=m['rexs'+str(s)+'-0']+m['rexs'+str(s)+'-2']+m['rexs'+str(s)+'-3']
        m['ims'+str(s)]=m['rims'+str(s)+'-0']+m['rims'+str(s)+'-2']+m['rims'+str(s)+'-3']
        m['trds'+str(s)]=m['exs'+str(s)]+m['ims'+str(s)]
        m['nxs'+str(s)]=m['nxs'+str(s)+'-0']+m['nxs'+str(s)+'-2']+m['nxs'+str(s)+'-3']

    m['ex_nafta']=m.rex0+m.rex2
    m['im_nafta']=m.rim0+m.rim2
    m['trd_nafta']=m.ex_nafta+m.im_nafta
    m['nx_nafta']=m.nx0+m.nx2
    m['reer_nafta']=(m.rer0*(m.rex0+m.rim0) + m.rer2*(m.rex2+m.rim2))/(m.rex0+m.rim0+m.rex2+m.rim2)

    for s in [0,1,2,3,4]:
        m['exs'+str(s)+'_nafta']=m['rexs'+str(s)+'-0']+m['rexs'+str(s)+'-2']
        m['ims'+str(s)+'_nafta']=m['rims'+str(s)+'-0']+m['rims'+str(s)+'-2']
        m['trds'+str(s)+'_nafta']=m['exs'+str(s)+'_nafta']+m['ims'+str(s)+'_nafta']
        m['nxs'+str(s)+'_nafta']=m['nxs'+str(s)+'-0']+m['nxs'+str(s)+'-2']
    
    m['trd_row']=m.rex3+m.rim3

for m in models['MEX']:
    m['nx']=m.nx0+m.nx1+m.nx3

    m['ex']=m.rex0+m.rex1+m.rex3
    m['im']=m.rim0+m.rim1+m.rim3
    m['trd']=m.ex+m.im
    m['reer']=(m.rer0*(m.rex0+m.rim0) + m.rer1*(m.rex1+m.rim1) + m.rer3*(m.rex3+m.rim3))/(m.rex0+m.rim0+m.rex1+m.rim1+m.rex3+m.rim3)

    for s in [0,1,2,3,4]:
        m['exs'+str(s)]=m['rexs'+str(s)+'-0']+m['rexs'+str(s)+'-1']+m['rexs'+str(s)+'-3']
        m['ims'+str(s)]=m['rims'+str(s)+'-0']+m['rims'+str(s)+'-1']+m['rims'+str(s)+'-3']
        m['trds'+str(s)]=m['exs'+str(s)]+m['ims'+str(s)]
        m['nxs'+str(s)]=m['nxs'+str(s)+'-0']+m['nxs'+str(s)+'-1']+m['nxs'+str(s)+'-3']

    m['ex_nafta']=m.rex0+m.rex1
    m['im_nafta']=m.rim0+m.rim1
    m['trd_nafta']=m.ex_nafta+m.im_nafta
    m['nx_nafta']=m.nx0+m.nx1
    m['reer_nafta']=(m.rer0*(m.rex0+m.rim0) + m.rer1*(m.rex1+m.rim1))/(m.rex0+m.rim0+m.rex1+m.rim1)

    for s in [0,1,2,3,4]:
        m['exs'+str(s)+'_nafta']=m['rexs'+str(s)+'-0']+m['rexs'+str(s)+'-1']
        m['ims'+str(s)+'_nafta']=m['rims'+str(s)+'-0']+m['rims'+str(s)+'-1']
        m['trds'+str(s)+'_nafta']=m['exs'+str(s)+'_nafta']+m['ims'+str(s)+'_nafta']
        m['nxs'+str(s)+'_nafta']=m['nxs'+str(s)+'-0']+m['nxs'+str(s)+'-1']

    m['trd_row']=m.rex3+m.rim3

################################################################################################

country_names = {'USA':'United States','CAN':'Canada','MEX':'Mexico','ROW':'rest of world','TOT':'world'}
country_names2 = {0:'United States',1:'Canada',2:'Mexico',3:'rest of world','TOT':'world'}
panels = {'USA':'(a)','CAN':'(b)','MEX':'(c)'}
countries=['USA','CAN','MEX']
partners={'USA':[1,2,3],'CAN':[0,2,3],'MEX':[0,1,3]}
partner_names={'USA':['CAN','MEX','ROW'],'CAN':['USA','MEX','ROW'],'MEX':['USA','CAN','ROW']}
sectors=[0,1,2,3,4]
NT=50
phi=1.0
TNAFTA=0


with open('output/model_lr_changes.tex','wb') as file:
    file.write('\\begin{table}[p]\n')
    file.write('\\renewcommand{\\arraystretch}{1.2}\n')
    file.write('\\begin{center}\n')
    file.write("\\caption{Long-run effects of NAFTA termination (percent changes)}\n")
    file.write('\\label{tab:lr_changes}\n')
    file.write('\\footnotesize\n')
    file.write('\\begin{tabular}{lcccccc}\n')
    
    file.write('\\toprule\n')
    file.write('\\multicolumn{1}{p{2cm}}{\\centering Quantity} & ')

    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Agriculture} & ')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Resource\\\\extraction} &')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Trans.\\\\equip.} &')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Other\\\\manuf.} & ')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Services} & ')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Total}\\\\\n')
    file.write('\\midrule\n')

    for c in countries:
        file.write('\\multicolumn{6}{l}{\\textit{'+panels[c]+' '+country_names[c]+'}}\\\\\n')

        file.write('Value added')
        for col in ['va0','va1','va2','va3','va4','rgdp']:
            val = 100.0*(models[c][1][col][NT]/models[c][0][col][NT]-1.0)
            file.write('& %0.2f' % val)
        file.write('\\\\\n')

        file.write('Consumption')
        for col in ['c0','c1','c2','c3','c4','c']:
            val = 100.0*(models[c][1][col][NT]/models[c][0][col][NT]-1.0)
            file.write('& %0.2f' % val)
        file.write('\\\\\n')

        file.write('Investment')
        for col in ['i0','i1','i2','i3','i4','ii']:
            val = 100.0*(models[c][1][col][NT]/models[c][0][col][NT]-1.0)
            file.write('& %0.2f' % val)
        file.write('\\\\\n')

        for p in ['']+partners[c]:
            if p=='':
                file.write('Exports')
            else:
                file.write('\quad to ' + country_names2[p])
            for s in ['0','1','2','3','4','']:
                col=''
                if s!='':
                    col='exs'+s
                    if p!='':
                        col=col+'-'+str(p)
                else:
                    col='ex'+str(p)
                  
                val = 100.0*(models[c][1][col][NT]/models[c][0][col][NT]-1.0)
                file.write('& %0.2f' % val)            
            file.write('\\\\\n')

        for p in ['']+partners[c]:
            if p=='':
                file.write('Imports')
            else:
                file.write('\quad from ' + country_names2[p])
            for s in ['0','1','2','3','4','']:
                col=''
                if s!='':
                    col='ims'+s
                    if p!='':
                        col=col+'-'+str(p)
                else:
                    col='im'+str(p)
                  
                val = 100.0*(models[c][1][col][NT]/models[c][0][col][NT]-1.0)
                file.write('& %0.2f' % val)            
            file.write('\\\\\n')

        if(c!='MEX'):
            file.write('\\\\\n')

    file.write('\\bottomrule\n')
    file.write('\\end{tabular}\n')
    file.write('\\normalsize\n')
    file.write('\\end{center}\n')
    file.write('\\end{table}\n')
