#####################################################################################
#imports, small functions, etc.

import numpy as np
from scipy import signal
import pandas as pd

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

phi=1.0
NT=50
TNAFTA=0

#################################################################################
# load the model results

def load_results(suff):

    models_usa = []
    models_usa.append(pd.read_csv('../c/output/vars0_usa'+suff+'.csv'))
    models_usa.append(pd.read_csv('../c/output/vars1_usa'+suff+'.csv'))
    for m in models_usa:
        m['im_nafta']=m.rim1+m.rim2
        m['ex_nafta']=m.rex1+m.rex2
        m['trd_nafta']=m.ex_nafta+m.im_nafta

    models_can = []
    models_can.append(pd.read_csv('../c/output/vars0_can'+suff+'.csv'))
    models_can.append(pd.read_csv('../c/output/vars1_can'+suff+'.csv'))
    for m in models_can:
        m['im_nafta']=m.rim0+m.rim2
        m['ex_nafta']=m.rex0+m.rex2
        m['trd_nafta']=m.ex_nafta+m.im_nafta


    models_mex = []
    models_mex.append(pd.read_csv('../c/output/vars0_mex'+suff+'.csv'))
    models_mex.append(pd.read_csv('../c/output/vars1_mex'+suff+'.csv'))
    for m in models_mex:
        m['im_nafta']=m.rim0+m.rim1
        m['ex_nafta']=m.rex0+m.rex1
        m['trd_nafta']=m.ex_nafta+m.im_nafta
        
    models={'USA':models_usa,'CAN':models_can,'MEX':models_mex}

    return models

models = load_results('')
models_s0 = load_results('_no_k_adj_cost')
models_s1 = load_results('_no_l_adj_cost')
models_s2 = load_results('_no_trd_adj_cost')
models_s3 = load_results('_eqkappa')
models_s4 = load_results('_nokappa')
models_s5 = load_results('_noio')
models_s6 = load_results('_sym_te')
models_s7 = load_results('_cobb_douglas')
models_s8 = load_results('_fix_tb')

models_a6 = load_results('_tariff_alt6')
models_s0a6 = load_results('_no_k_adj_cost_iceberg')
models_s1a6 = load_results('_no_l_adj_cost_iceberg')
models_s2a6 = load_results('_no_trd_adj_cost_iceberg')
models_s3a6 = load_results('_eqkappa_iceberg')
models_s4a6 = load_results('_nokappa_iceberg')
models_s5a6 = load_results('_noio_iceberg')
models_s6a6 = load_results('_sym_te_iceberg')
models_s7a6 = load_results('_cobb_douglas_iceberg')

models_a1 = load_results('_tariff_alt1')
models_a4 = load_results('_tariff_alt4')
models_a4_v2 = load_results('_tariff_alt4_v2')
models_a5 = load_results('_tariff_alt5')
models_a2 = load_results('_tariff_alt2')
models_a7 = load_results('_tariff_alt7')
models_cp = load_results('_cp_combo')


#################################################################################
# write the table

file=open('output/welfare.tex','wb')

#file.write('\\begin{landscape}\n')
file.write('\\begin{table}[h!]\n')
file.write('\\footnotesize\n')
file.write('\\renewcommand{\\arraystretch}{1.2}\n')
file.write('\\begin{center}\n')

file.write('\\caption{Welfare effects of NAFTA termination}\n')
file.write('\\label{tab:welfare}\n')

# number of columns
file.write('\\begin{tabular}{lcccccccccc}')
file.write('\\toprule\n')

# headers
file.write('& \\multicolumn{3}{c}{Dynamic (pct. change)} & \\multicolumn{3}{c}{Long-run (pct. change)} & \\multicolumn{3}{c}{Ratio dynamic to long-run}\\\\\n')
file.write('\\cmidrule(rl){2-4}\\cmidrule(rl){5-7}\\cmidrule(rl){8-10}\n')
file.write('Model')
for c in ['USA','Canada','Mexico']:
    file.write('& \\multicolumn{1}{p{1.0cm}}{\\centering '+c+'}')
for c in ['USA','Canada','Mexico']:
    file.write('& \\multicolumn{1}{p{1.0cm}}{\\centering '+c+'}')
for c in ['USA','Canada','Mexico']:
    file.write('& \\multicolumn{1}{p{1.0cm}}{\\centering '+c+'}')
file.write('\\\\\n\\midrule\n')

def write_results(mm,label):
    file.write(label)
    for c in ['USA','CAN','MEX']:
        tmp = 100.0*((mm[c][1].W[TNAFTA]/mm[c][0].W[TNAFTA])**(1.0/phi)-1.0)
        file.write('& %0.3f' % tmp)
    for c in ['USA','CAN','MEX']:
        tmp = 100.0*(mm[c][1].c[NT]/mm[c][0].c[NT]-1.0)
        file.write('& %0.3f' % tmp)
    for c in ['USA','CAN','MEX']:
        tmp1 = 100.0*((mm[c][1].W[TNAFTA]/mm[c][0].W[TNAFTA])**(1.0/phi)-1.0)
        tmp2 = 100.0*(mm[c][1].c[NT]/mm[c][0].c[NT]-1.0)
        tmp3 = tmp1/tmp2
        file.write('& %0.3f' % tmp3)
    file.write('\\\\\n')

write_results(models,'Baseline')

file.write('\\\\\n')
file.write('\\multicolumn{10}{l}{\\textit{(a) Effects of dynamic ingredients}}')
file.write('\\\\\n')

write_results(models_s0,'No capital adj. costs')
write_results(models_s1,'No labor adj. costs')
write_results(models_s2,'No import adj. costs')
write_results(models_s3,'Static exporting')
write_results(models_s4,'No extensive margin')
write_results(models_s8,'Fixed trade balances')

file.write('\\\\\n')
file.write('\\multicolumn{10}{l}{\\textit{(b) Effects of static ingredients}}')
file.write('\\\\\n')
write_results(models_s5,'No intermediate inputs')
write_results(models_s7,'Cobb-Douglas production')
write_results(models_s6,'Sym. trade elasticities')

# file.write('\\\\\n')
# file.write('\\multicolumn{10}{l}{\\textit{Iceberg trade costs and\ldots}}')
# file.write('\\\\\n')
# write_results(models_a6,'Baseline model')
# write_results(models_s0a6,'No capital adj. costs')
# write_results(models_s1a6,'No labor adj. costs')
# write_results(models_s2a6,'No import adj. costs')
# write_results(models_s3a6,'Static exporting')
# write_results(models_s4a6,'No extensive margin')
# write_results(models_s5a6,'No intermediate inputs')
# write_results(models_s6a6,'Sym. trade elasticities')
# write_results(models_s7a6,'Cobb-Douglas production')

file.write('\\\\\n')
file.write('\\multicolumn{10}{l}{\\textit{(c) Alternative scenarios}}')
file.write('\\\\\n')

write_results(models_a4,'USMCA')
write_results(models_a4_v2,'Stricter dom. content reqs.')
write_results(models_a5,'US-Canada FTA')
write_results(models_a2,'Canada-Mexico FTA')
write_results(models_a1,'Higher U.S. tariffs')
write_results(models_a7,'1994 tariffs')
write_results(models_cp,'CP specification')
#write_results(models_a6,'Iceberg costs, not tariffs')

file.write('\\bottomrule\n')
file.write('\\end{tabular}\n')
file.write('\\end{center}\n')
file.write('\\normalsize\n')
file.write('\\end{table}\n')
#file.write('\\end{landscape}\n')
file.close()

