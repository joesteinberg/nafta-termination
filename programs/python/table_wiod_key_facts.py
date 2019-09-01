import numpy as np
import pandas as pd

##################################################################################
# setup
m = pd.read_pickle('output/wiod_m.pik')
f = pd.read_pickle('output/wiod_f.pik')
vg = pd.read_pickle('output/wiod_vg.pik')

##################################################################################
# calculate trade balance and investment rate data

# intermediate trade
m_trd =  m.groupby(['year','col_region','row_region','row_sector'])['M'].sum().reset_index()
m_trd=m_trd[m_trd.col_region != m_trd.row_region]
im_m = m_trd.rename(columns={'col_region':'region','row_sector':'sector','row_region':'partner','M':'im_M'})
ex_m = m_trd.rename(columns={'row_region':'region','row_sector':'sector','col_region':'partner','M':'ex_M'})
m_trd2 = pd.merge(left=ex_m,right=im_m,how='left',on=['year','region','partner','sector'])

# final trade
f_trd = f[f.col_region != f.row_region]
im_f = f_trd.rename(columns={'col_region':'region','row_sector':'sector','row_region':'partner','C':'im_C','I':'im_I'})
ex_f = f_trd.rename(columns={'row_region':'region','row_sector':'sector','col_region':'partner','C':'ex_C','I':'ex_I'})
f_trd2 = pd.merge(left=ex_f,right=im_f,how='left',on=['year','region','partner','sector'])

# merge and calculate totals + balances
trd = pd.merge(left=m_trd2,right=f_trd2,how='left',on=['year','region','partner','sector'])

for d in ['ex','im']:
    trd[d+'_F'] = trd[d+'_C']+trd[d+'_I']
    trd[d] = trd[d+'_M']+trd[d+'_F']

for u  in ['_M','_C','_I','_F','']:
    trd['tb'+u] = trd['ex'+u] - trd['im'+u]

# aggregate by sector and append
cols = []
for d in ['ex','im','tb']:
    for u in ['','_M','_F','_C','_I']:
        cols.append(d+u)

g = trd.groupby(['year','region','partner'])
sums = g[cols].sum().reset_index()
sums['sector'] = 'TOT'
trd = trd.append(sums)

# aggregate by country and append
g = trd.groupby(['year','region','sector'])
sums = g[cols].sum().reset_index()
sums['partner'] = 'TOT'
trd = trd.append(sums)
trd = trd.sort_values(['year','region','partner','sector']).reset_index(drop=True)

# merge on value added
va = vg.groupby(['year','col_region','col_sector'])['VA'].sum().reset_index()
va.rename(columns={'col_region':'region','col_sector':'sector'},inplace=True)
trd = pd.merge(left=trd,right=va,how='left',on=['year','region','sector'])

# merge on consumption
cons = f.groupby(['year','col_region','row_sector'])['C'].sum().reset_index()
cons.rename(columns={'col_region':'region','row_sector':'sector'},inplace=True)
trd = pd.merge(left=trd,right=cons,how='left',on=['year','region','sector'])

# merge on gdp
gdp = vg.groupby(['year','col_region'])['VA'].sum().reset_index()
gdp.rename(columns={'col_region':'region','VA':'GDP'},inplace=True)
trd = pd.merge(left=trd,right=gdp,how='left',on=['year','region'])
trd.loc[trd.sector=='TOT','VA']=trd.loc[trd.sector=='TOT','GDP']

# save
trd.to_pickle('output/wiod_trd.pik')

# make latex table summarizing 2014 data
#trd=trd[trd.year.isin(range(2010,2015))].groupby(['region','partner','sector']).mean().reset_index()
trd=trd[trd.year==2014].groupby(['region','partner','sector']).mean().reset_index()

sector_names = {'A':'Agriculture',
                'R':'Resource extraction',
                'T':'Transport equipment',
                'M':'Other manufacturing',
                'S':'Services',
                'TOT':'All goods'}

country_names = {'USA':'United States','CAN':'Canada','MEX':'Mexico','ROW':'rest of world','TOT':'world'}
partners = {'USA':['TOT','CAN','MEX','ROW'],'CAN':['TOT','USA','MEX','ROW'],'MEX':['TOT','USA','CAN','ROW']}
panels = {'USA':'(a)','CAN':'(b)','MEX':'(c)'}

with open('output/key_facts.tex','wb') as file:
    file.write('\\begin{table}[p]\n')
    file.write('\\renewcommand{\\arraystretch}{1.2}\n')
    file.write('\\begin{center}\n')
    file.write("\\caption{Sectoral production and trade in NAFTA (2014 data, percent GDP)}\n")
    file.write('\\label{tab:key_facts}\n')
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

    for c in ['USA','CAN','MEX']:
        file.write('\\multicolumn{6}{l}{\\textit{'+panels[c]+' '+country_names[c]+'}}\\\\\n')
        mask=trd.region==c

        file.write('Value added')
        mask2 = np.logical_and(mask,trd.partner=='TOT')
        for s in ['A','R','T','M','S','TOT']:
            mask3=np.logical_and(mask2,trd.sector==s)
            masked=trd[mask3]
            val = 100.0*masked['VA']/masked['GDP']
            file.write('& %0.2f' % val)
        file.write('\\\\\n')
        
        for p in partners[c]:
            mask2=np.logical_and(mask,trd.partner==p)
            if p=='TOT':
                file.write('Exports')
            else:
                file.write('\quad to ' + country_names[p])
            for s in ['A','R','T','M','S','TOT']:
                mask3=np.logical_and(mask2,trd.sector==s)
                masked=trd[mask3]
                val = 100.0*(masked['ex'])/masked['GDP']
                file.write('& %0.2f' % val)            
            file.write('\\\\\n')

        for p in partners[c]:
            mask2=np.logical_and(mask,trd.partner==p)
            if p=='TOT':
                file.write('Imports')
            else:
                file.write('\quad from ' + country_names[p])
            for s in ['A','R','T','M','S','TOT']:
                mask3=np.logical_and(mask2,trd.sector==s)
                masked=trd[mask3]
                val = 100.0*(masked['im'])/masked['GDP']
                file.write('& %0.2f' % val)            
            file.write('\\\\\n')

        for p in partners[c]:
            mask2=np.logical_and(mask,trd.partner==p)
            if p=='TOT':
                file.write('Net exports')
            else:
                file.write('\quad with ' + country_names[p])
            for s in ['A','R','T','M','S','TOT']:
                mask3=np.logical_and(mask2,trd.sector==s)
                masked=trd[mask3]
                val = 100.0*masked['tb']/masked['GDP']
                file.write('& %0.2f' % val)            
            file.write('\\\\\n')

        if(c!='MEX'):
            file.write('\\\\\n')

    file.write('\\bottomrule\n')
    file.write('\\end{tabular}\n')
    file.write('\\normalsize\n')
    file.write('\\end{center}\n')
    file.write('\\end{table}\n')


