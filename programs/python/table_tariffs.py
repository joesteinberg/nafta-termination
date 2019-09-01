#############################################################################
# This script calculates MFN tariff rates in bilateral trade between NAFTA
# members. MFN tariffs are reported at the 6-digit HS level for each country.
# I merge these with bilateral trade flows from COMTRADE to use as weights.
# I assign each HS industry to one of the 4 goods sectors in the input-output
# table (agriculture, resource extraction, transportation equipment, other
# manufacturing); services has no tariffs.
#############################################################################

import numpy as np
import pandas as pd

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
        return ''

d_country_code = {'USA':0,'CAN':1,'MEX':2,'WLD':3}
d_sector_code = {'A':0,'R':1,'T':2,'M':3,'TOT':4}

sector_names = {'A':'Agriculture',
                'R':'Resource extraction',
                'T':'Transport equipment',
                'M':'Other manufacturing',
                'TOT':'All goods'}

country_names = {'USA':'United States','CAN':'Canada','MEX':'Mexico','WLD':'Rest of world'}
panels = {'USA':'(a)','CAN':'(b)','MEX':'(c)'}

def add0(hs12):
    x=len(hs12)
    if x==5:
        return '0'+hs12
    else:
        return hs12

def country_code(country):
    if country=='Canada':
        return 'CAN'
    elif (country=='United States of America' or country=='United States'):
        return 'USA'
    elif country=='Mexico':
        return 'MEX'
    else:
        return ''


##################################################################3
# load tariff data

trains = pd.read_csv('../../data/trains_data_extract.csv')
trains = trains.replace('..',np.nan)
trains = trains[~(trains['Product Code'].isnull())]
trains['country']=trains.Reporter.apply(country_code)
trains['partner']=trains['Partner Name'].apply(country_code)
trains['hs_code']=trains['Product Code'].apply(lambda x: add0(str(int(x))))
trains = trains[trains['Tariff Indicator Code']=='SimpleAverage']

trains_mex = trains[['country','partner','hs_code','1991 [1991]']][trains.country=='MEX']
trains_mex = trains_mex.rename(columns={'1991 [1991]':'avg_duty'})
trains_mex = trains_mex[~(trains_mex.avg_duty.isnull())]

trains_other = trains[['country','partner','hs_code','1993 [1993]']][trains.country.isin(['USA','CAN'])]
trains_other = trains_other.rename(columns={'1993 [1993]':'avg_duty'})
trains_other = trains_other[~(trains_other.avg_duty.isnull())]

trains2 = trains_other.append(trains_mex)

##################################################################3
# process 2014 data

comtrade = pd.read_csv('../../data/comtrade_nafta_2014_hs12.csv',
                       dtype={'country':str,
                              'partner':str,
                              'hs_code':str,
                              'imports':float})

comtrade['hs_code']=comtrade.hs_code.apply(add0)
comtrade = comtrade[comtrade.hs_code != '999999']
comtrade = comtrade[comtrade.country != comtrade.partner]

sums_nafta = comtrade[comtrade.partner!='WLD'].groupby(['country','hs_code'])['imports'].sum()
sums_nafta=sums_nafta.reset_index().rename(columns={'imports':'nafta_imports'})
comtrade = pd.merge(left=comtrade,right=sums_nafta,how='left',on=['country','hs_code'])
comtrade.nafta_imports.fillna(0,inplace=True)
comtrade.loc[comtrade.partner=='WLD','imports'] = comtrade.loc[comtrade.partner=='WLD','imports']-comtrade.loc[comtrade.partner=='WLD','nafta_imports']

wto = pd.read_csv('../../data/wto_mfn_tariffs_nafta_2014_hs12.csv',
                       dtype={'country':str,
                              'hs_level':int,
                              'hs_code':str,
                              'avg_duty':float})


wto['country'] = wto.country.apply(country_code)
wto = wto[wto.hs_level==6]
wto.avg_duty[wto.avg_duty.isnull()] = 0.0
wto=wto.drop('hs_level',axis=1)


merged = pd.merge(left=comtrade,
                  right=wto,
                  how='left',
                  on=['country','hs_code'])

merged['hs2'] = merged.hs_code.apply(lambda x: hs2(x))
merged['sector'] = merged.hs2.apply(lambda x: sector(x))

g = merged.groupby(['country','partner','sector'])

wavg = lambda x: np.average(x,weights=merged.loc[x.index,'imports'])
result  = g['avg_duty'].agg(wavg).reset_index()

result['country_code']=result.country.apply(lambda x: d_country_code[x])
result['partner_code']=result.partner.apply(lambda x: d_country_code[x])
result['sector_code']=result.sector.apply(lambda x: d_sector_code[x])
result.sort_values(by=['country_code','partner_code','sector_code'],inplace=True)
result=result.reset_index(drop=True)

# write the results to txt and latex... in both cases, we set the tariff on ROW trade to zero
# we want all tariffs to be zero in period 0 to faciliate calibration (See paper text)
result_2014_csv = result[['country_code','partner_code','sector_code','avg_duty']]
result_2014_csv.loc[result_2014_csv.partner_code==3,'avg_duty']=0.0
result_2014_csv.to_csv('output/tariffs.txt',
              columns=['country_code','partner_code','sector_code','avg_duty'],
              index=False,
              header=False,
              sep=' ')

g2 = merged.groupby(['country','partner'])
result2 = g2['avg_duty'].agg(wavg).reset_index()
result2['sector']='TOT'
result2['country_code']=result2.country.apply(lambda x: d_country_code[x])
result2['partner_code']=result2.partner.apply(lambda x: d_country_code[x])
result2['sector_code']=result2.sector.apply(lambda x: d_sector_code[x])
result2.sort_values(by=['country_code','partner_code','sector_code'],inplace=True)
result2=result2.reset_index(drop=True)

result_2014_latex=result.append(result2)
result_2014_latex.sort_values(by=['country_code','partner_code','sector_code'],inplace=True)
result_2014_latex=result_2014_latex.reset_index(drop=True)
result_2014_latex.loc[result_2014_latex.partner=='WLD','avg_duty']=0.0

with open('output/tariffs.tex','wb') as file:
    file.write('\\begin{table}[p]\n')
    file.write('\\renewcommand{\\arraystretch}{1.2}\n')
    file.write('\\begin{center}\n')
    file.write('\\caption{Change in import tariffs after NAFTA termination}\n')
    file.write('\\label{tab:mfn_tariffs}\n')
    file.write('\\footnotesize\n')
    file.write('\\begin{tabular}{lccccc}\n')    
    file.write('\\toprule\n')
    #file.write('\\multicolumn{1}{p{2cm}}{\\centering Country} &')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Trade\\\\partner} & ')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Agriculture} & ')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Resource\\\\extraction} &')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Trans.\\\\equip.} &')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Other\\\\manuf.} & ')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Total}\\\\\n')
    file.write('\\midrule\n')

    for c in ['USA','CAN','MEX']:
        file.write('\\multicolumn{6}{l}{\\textit{'+panels[c]+' '+country_names[c]+'}}\\\\\n')
        mask=result_2014_latex.country==c
        for p in ['USA','CAN','MEX']:
            if(p!=c):
                mask2=np.logical_and(mask,result_2014_latex.partner==p)
                #file.write(c + '&' + p)
                file.write('\\quad '+country_names[p])
                for s in ['A','R','T','M','TOT']:
                    mask3=np.logical_and(mask2,result_2014_latex.sector==s)
                    x= result_2014_latex.loc[mask3,'avg_duty'].values[0]
                    file.write('& %0.2f' % x)
                file.write('\\\\\n')
        if(c!='MEX'):
            file.write('\\\\\n')

    file.write('\\bottomrule\n')
    file.write('\\end{tabular}\n')
    file.write('\\normalsize\n')
    file.write('\\end{center}\n')
    file.write('\\end{table}\n')

    
##################################################################3
# alternative with higher US tariffs

result_alt1_csv = result.copy()
result_alt1_csv.loc[np.logical_and(result_alt1_csv.country=='USA',result_alt1_csv.partner != 'WLD'),'avg_duty'] = result_alt1_csv.loc[np.logical_and(result_alt1_csv.country=='USA',result_alt1_csv.partner != 'WLD'),'avg_duty']*2.0
result_alt1_csv.loc[np.logical_and(result_alt1_csv.country!='USA',result_alt1_csv.partner == 'WLD'),'avg_duty'] = 0.0

result_alt1_csv.to_csv('output/tariffs_alt1.txt',
               columns=['country_code','partner_code','sector_code','avg_duty'],
               index=False,
               header=False,
               sep=' ')

##################################################################3
# process 1993 data

comtrade_old = pd.read_csv('../../data/comtrade_nafta_2014_hs96.csv',
                           dtype={'country':str,
                                  'partner':str,
                                  'hs_code':str,
                                  'imports':float})

comtrade_old['hs_code']=comtrade_old.hs_code.apply(add0)
comtrade_old = comtrade_old[comtrade_old.hs_code != '999999']
comtrade_old = comtrade_old[comtrade_old.country != comtrade_old.partner]

tmp = comtrade_old[comtrade_old.imports.isnull()]['hs_code'].unique()
tmp2 = comtrade[['country','partner','hs_code','imports']][comtrade.hs_code.isin(tmp)].rename(columns={'imports':'imports_'})
comtrade_old = pd.merge(left=comtrade_old,right=tmp2,how='left',on=['country','partner','hs_code'])
comtrade_old.loc[comtrade_old.imports.isnull(),'imports'] = comtrade_old.loc[comtrade_old.imports.isnull(),'imports_']
comtrade_old.drop('imports_',axis=1)


sums_nafta = comtrade_old[comtrade_old.partner!='WLD'].groupby(['country','hs_code'])['imports'].sum()
sums_nafta=sums_nafta.reset_index().rename(columns={'imports':'nafta_imports'})
comtrade_old = pd.merge(left=comtrade_old,right=sums_nafta,how='left',on=['country','hs_code'])
comtrade_old.nafta_imports.fillna(0,inplace=True)
comtrade_old.loc[comtrade_old.partner=='WLD','imports'] = comtrade_old.loc[comtrade_old.partner=='WLD','imports']-comtrade_old.loc[comtrade_old.partner=='WLD','nafta_imports']

merged = pd.merge(left=comtrade_old,
                  right=trains2,
                  how='left',
                  on=['country','partner','hs_code'])
merged=merged[~(merged.avg_duty.isnull())].reset_index()
merged['avg_duty'] = merged.avg_duty.astype(float)
merged['hs2'] = merged.hs_code.apply(lambda x: hs2(x))
merged['sector'] = merged.hs2.apply(lambda x: sector(x))

g = merged.groupby(['country','partner','sector'])

wavg = lambda x: np.average(x,weights=merged.loc[x.index,'imports'])
result  = g['avg_duty'].agg(wavg).reset_index()

result['country_code']=result.country.apply(lambda x: d_country_code[x])
result['partner_code']=result.partner.apply(lambda x: d_country_code[x])
result['sector_code']=result.sector.apply(lambda x: d_sector_code[x])
result.sort_values(by=['country_code','partner_code','sector_code'],inplace=True)
result=result.reset_index(drop=True)


result_1993_csv = result[['country_code','partner_code','sector_code','avg_duty']]
for i in range(3):
    for s in range(4):
        result_1993_csv = result_1993_csv.append({'country_code':i,
                                                  'partner_code':3,
                                                  'sector_code':s,
                                                  'avg_duty':0},ignore_index=True)
result_1993_csv = result_1993_csv.sort_values(by=['country_code','partner_code','sector_code'],
                                              ascending=[True,True,True])

result_1993_csv.to_csv('output/tariffs_old.txt',
              columns=['country_code','partner_code','sector_code','avg_duty'],
              index=False,
              header=False,
              sep=' ')

g2 = merged.groupby(['country','partner'])
result2 = g2['avg_duty'].agg(wavg).reset_index()
result2['sector']='TOT'
result2['country_code']=result2.country.apply(lambda x: d_country_code[x])
result2['partner_code']=result2.partner.apply(lambda x: d_country_code[x])
result2['sector_code']=result2.sector.apply(lambda x: d_sector_code[x])
result2.sort_values(by=['country_code','partner_code','sector_code'],inplace=True)
result2=result2.reset_index(drop=True)

result_1993_latex=result.append(result2)
result_1993_latex.sort_values(by=['country_code','partner_code','sector_code'],inplace=True)
result_1993_latex=result_1993_latex.reset_index(drop=True)
result_1993_latex.loc[result_1993_latex.partner=='WLD','avg_duty']=0.0

with open('output/tariffs_old.tex','wb') as file:
    file.write('\\begin{table}[p]\n')
    file.write('\\renewcommand{\\arraystretch}{1.2}\n')
    file.write('\\begin{center}\n')
    file.write('\\caption{Change in import tariffs to pre-NAFTA levels}\n')
    file.write('\\label{tab:mfn_tariffs_old}\n')
    file.write('\\footnotesize\n')
    file.write('\\begin{tabular}{lccccc}\n')    
    file.write('\\toprule\n')
    #file.write('\\multicolumn{1}{p{2cm}}{\\centering Country} &')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Trade\\\\partner} & ')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Agriculture} & ')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Resource\\\\extraction} &')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Trans.\\\\equip.} &')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Other\\\\manuf.} & ')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Total}\\\\\n')
    file.write('\\midrule\n')

    for c in ['USA','CAN','MEX']:
        file.write('\\multicolumn{6}{l}{\\textit{'+panels[c]+' '+country_names[c]+'}}\\\\\n')
        mask=result_1993_latex.country==c
        for p in ['USA','CAN','MEX']:
            if(p!=c):
                mask2=np.logical_and(mask,result_1993_latex.partner==p)
                #file.write(c + '&' + p)
                file.write('\\quad '+country_names[p])
                for s in ['A','R','T','M','TOT']:
                    mask3=np.logical_and(mask2,result_1993_latex.sector==s)
                    x= result_1993_latex.loc[mask3,'avg_duty'].values[0]
                    file.write('& %0.2f' % x)
                file.write('\\\\\n')
        if(c!='MEX'):
            file.write('\\\\\n')

    file.write('\\bottomrule\n')
    file.write('\\end{tabular}\n')
    file.write('\\normalsize\n')
    file.write('\\end{center}\n')
    file.write('\\end{table}\n')
