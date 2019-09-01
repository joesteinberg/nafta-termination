##########################################################################################
# This script calculates trade elasticities for US, CAN, MEX, and ROW at the goods-sector level,
# for intermediate and final trade separately. Results are averages of Caliendo and Parro (2014)
# elasctities for 2-digit ISIC sectors (the same sectors as in the WIOD data) weighted by
# imports in the 2014 WIOD data.
##########################################################################################

import numpy as np
import pandas as pd

year = 2014

# country list (aside from ROW, which we will construct by aggregating all other countries)
regions = ['USA','CAN','MEX']
region_num = {'USA':0,'CAN':1,'MEX':2,'ROW':3}

#########################################################################################
# baseline CP (2015) elasticities

# sector aggregation scheme
# A: agriculture
# R: resource extraction
# T: cars + car parts, other transportation equipment
# M: other mfg
# S: services
sectors={'A':[1,2,3],'R':[4,10],'T':[20,21],'M':range(5,10)+range(11,20)+[22],'S':range(23,57)}
sector_num={'A':0,'R':1,'T':2,'M':3,'S':4}

# elasticities for industries from Caliendo and Parro
# notes:
# 1. all 3 agriculture sectors (crop and animals, forestry and logging, fishing) get CP's ag value of 8.11
# 2. both paper and printing get CP's paper+printing value of 9.07
# 3. both chemicals and pharma get CP's chem value of 4.75
# 4. computer, electronic, and optical products get average of CP's values for office, communication, and medical)
elasticities = {1:8.11,
                2:8.11,
                3:8.11,
                4:15.72,
                5:2.55,
                6:5.56,
                7:10.83,
                8:9.07,
                9:9.07,
                10:51.08,
                11:4.75,
                12:4.75,
                13:1.66,
                14:2.76,
                15:7.99,
                16:4.30,
                17:(12.79+7.07+8.98)/3.0,
                18:10.60,
                19:1.52,
                20:1.01,
                21:0.37,
                22:5.0}

def assign_elasticity(industry):
    if industry <= 22:
        return elasticities[industry]
    else:
        return 5.0

# set sector and region codes
def which_sector(x):
    for i in sectors.keys():
        if x in sectors[i]:
            return i
    return -1

def which_region(c):
    if c == 'TOT':
        return 'TOT'
    else:
        for i in regions:
            if c==i:
                return i
        return 'ROW'

def which_use(u):
    if u in range(1,56):
        return 'M'
    elif u in [57,58,59]:
        return 'C'
    elif u in [60]:
        return 'I'

w = pd.read_pickle('/home/joe/Datasets/WIOD2016/stata/wiot_stata_Nov16/wiot_2014.pik')
w['row_sector']=w['row_code'].apply(which_sector)
w['col_sector']=w['col_code'].apply(which_sector)
w['row_region']=w['row_country'].apply(which_region)
w['col_region']=w['col_country'].apply(which_region)
w['col_use']=w['col_code'].apply(which_use)
w['row_e']=w['row_code'].apply(assign_elasticity)

wavg = lambda x: np.average(x,weights=w.loc[x.index,'value'])

# intermediate elasticities
#mask = w.col_use=='M'
#mask = np.logical_and(mask,w.col_region.isin(regions+['ROW'])) # ... valid use region
#mask = np.logical_and(mask,w.row_sector>=0)  # ...valid source sector
#mask = np.logical_and(mask,w.row_region.isin(regions+['ROW'])) # ... valid source region
## group by and aggregate
#g = w[mask].groupby(['col_region','row_sector'])
#m = g['row_e'].aggregate(wavg).reset_index()
#m.rename(columns={'col_region':'region','row_sector':'sector','row_e':'intermediate_elasticity'},inplace=True)

# aggregate elasticities
mask = w.col_use.isin(['M','C','I'])
mask = np.logical_and(mask,w.col_region.isin(regions+['ROW'])) # ... valid use region
mask = np.logical_and(mask,w.row_sector>=0)  # ...valid source sector
mask = np.logical_and(mask,w.row_region.isin(regions+['ROW'])) # ... valid source region
# group by and aggregate
g = w[mask].groupby(['col_region','row_sector'])
e = g['row_e'].aggregate(wavg).reset_index()
e.rename(columns={'col_region':'region','row_sector':'sector','row_e':'elasticity'},inplace=True)
e2 = e.groupby('sector')['elasticity'].mean().reset_index()

#e = pd.merge(left=m,right=f,how='left',on=['region','sector'])
#e['region_num'] = e.region.apply(lambda x: region_num[x])
#e['sector_num'] = e.sector.apply(lambda x: sector_num[x])
#e = e.sort_values(by=['region_num','sector_num'])
#e.to_csv('output/elasticities.txt',
#         columns=['elasticity'],
#         index=False,
#         header=False,
#         sep=' ')

#country_names = {'USA':'United States','CAN':'Canada','MEX':'Mexico','ROW':'Rest of world'}
#panels = {'USA':'(a)','CAN':'(b)','MEX':'(c)','ROW':'(d)'}

# with open('output/elasticities.tex','wb') as file:
#     file.write('\\begin{table}[p]\n')
#     file.write('\\renewcommand{\\arraystretch}{1.2}\n')
#     file.write('\\begin{center}\n')
#     file.write("\\caption{Assigned long-run trade elasticities}\n")
#     file.write('\\label{tab:elasticities}\n')
#     file.write('\\footnotesize\n')
#     file.write('\\begin{tabular}{lccccc}\n')
    
#     file.write('\\toprule\n')
#     file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Country} & ')
#     file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Agriculture} & ')
#     file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Resource\\\\extraction} &')
#     file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Trans.\\\\equip.} &')
#     file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Other\\\\manuf.} & ')
#     file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Services}\\\\\n')
#     file.write('\\midrule\n')

#     for c in ['USA','CAN','MEX','ROW']:
#         #file.write('\\multicolumn{6}{l}{\\textit{'+panels[c]+' '+country_names[c]+'}}\\\\\n')
#         mask=e.region==c

# #        file.write('Intermediate')
#         for s in ['A','R','T','M','S']:
#             mask2=np.logical_and(mask,e.sector==s)
#             x= e.loc[mask2,'intermediate_elasticity'].values[0]
#             file.write('& %0.2f' % x)
#         file.write('\\\\\n')

#         file.write('Final')
#         for s in ['A','R','T','M','S']:
#             mask2=np.logical_and(mask,e.sector==s)
#             x= e.loc[mask2,'final_elasticity'].values[0]
#             file.write('& %0.2f' % x)
#         file.write('\\\\\n')
        
#         if(c!='ROW'):
#             file.write('\\\\\n')

#     file.write('\\bottomrule\n')
#     file.write('\\end{tabular}\n')
#     file.write('\\normalsize\n')
#     file.write('\\end{center}\n')
#     file.write('\\end{table}\n')

##############################################################################################
# LTP version

def hs2(hs12):
    x=len(hs12)
    y=0
    if x==5:
        y=int(hs12[0])
    else:
        y=int(hs12[0:2])
    return y
    
def sector(hs2):

    if hs2<=14:
        return 'A'
    elif hs2<=24:
        return 'M'
    elif hs2<=27:
        return 'R'
    elif hs2<=38:
        return 'M'
    elif hs2<=40:
        return 'M'
    elif hs2<=43:
        return 'M'
    elif hs2<=46:
        return 'M'
    elif hs2<=49:
        return 'M'
    elif hs2<=63:
        return 'M'
    elif hs2<=67:
        return 'M'
    elif hs2<=70:
        return 'M'
    elif hs2==71:
        return 'M'
    elif hs2<=83:
        return 'M'
    elif hs2<=85:
        return 'M'
    elif hs2<=89:
        return 'T'
    elif hs2<=92:
        return 'M'
    elif hs2==93:
        return 'M'
    elif hs2<=96:
        return 'M'
    elif hs2==97:
        return 'M'
    else:
        hs2=''

def add0(hs12):
    x=len(hs12)
    if x==5:
        return '0'+hs12
    else:
        return hs12

def country_code(country):
    if country=='Canada':
        return 'CAN'
    elif country=='United States of America':
        return 'USA'
    elif country=='Mexico':
        return 'MEX'
    else:
        return ''

comtrade = pd.read_csv('../../data/comtrade_nafta_2014_hs12.csv',
                       dtype={'country':str,
                              'partner':str,
                              'hs_code':str,
                              'imports':float})

comtrade = comtrade[comtrade.partner=='WLD']
comtrade.drop('partner',axis=1,inplace=True)

comtrade['hs_code']=comtrade.hs_code.apply(add0)
comtrade = comtrade[comtrade.hs_code != '999999']

hs =comtrade.hs_code.unique().tolist()
countries=comtrade.country.unique().tolist()
index=pd.MultiIndex.from_product([countries,hs],names=['country','hs_code'])
combos=pd.DataFrame(index=index).reset_index()
comtrade = pd.merge(left=combos,right=comtrade,how='left',on=['country','hs_code'])
comtrade['imports']=comtrade.imports.fillna(0.0)

comtrade['hs2'] = comtrade.hs_code.apply(lambda x: hs2(x))
comtrade['sector'] = comtrade.hs2.apply(lambda x: sector(x))
totals = comtrade.groupby(['country'])['imports'].sum().reset_index().rename(columns={'imports':'total'})
comtrade = pd.merge(left=comtrade,right=totals,how='left',on=['country'])
comtrade['share']=comtrade.imports/comtrade.total
comtrade.drop('total',axis=1,inplace=True)

comtrade.sort_values(by=['country','share'],ascending=[True,True],inplace=True)
comtrade.reset_index(inplace=True,drop=True)
comtrade['cumshare']=comtrade.groupby(['country'])['share'].transform(np.cumsum)
comtrade['decile']=comtrade.cumshare.apply(lambda s: min((int(s*10))+1,10))

stotals = comtrade.groupby(['country','sector'])['imports'].sum().reset_index().rename(columns={'imports':'stotal'})
comtrade = pd.merge(left=comtrade,right=stotals,how='left',on=['country','sector'])
comtrade['sshare']=comtrade.imports/comtrade.stotal
comtrade.drop('stotal',axis=1,inplace=True)
ltp_share = comtrade[comtrade.decile==1].groupby(['country','sector'])['sshare'].sum().reset_index()
ls2= ltp_share.groupby('sector')['sshare'].mean().reset_index()

merged = pd.merge(left='e2',right='ls2',how='left',on='sector')
merged['e_ltp']=merged.e*(1.0+3.65*merged.ls2)

print merged




    



