import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

##################################################################################
# setup
print 'Performing setup tasks...'

import pandas as pd

years = [2014]

for y in years:
    stata = pd.read_stata('../data/WIOT'+str(y)+'_October16_ROW.dta')
    stata.drop(['IndustryCode','IndustryDescription','Year'],axis=1,inplace=True)
    stata.rename(columns={'Country':'row_country','RNr':'row_code'},inplace=True)
    cols = stata.columns.tolist()
    cols = cols[0:2]+[c[1:] for c in cols[2:-1]] + [cols[-1]+'62']
    stata.columns=cols
    stata2 = pd.melt(stata,id_vars=cols[0:2],value_vars=cols[2:])
    stata2['col_country'] = stata2.variable.str.slice(0,3,1)
    stata2['col_code'] = stata2.variable.str.slice(3,None,1).astype(int)
    stata2.drop('variable',axis=1,inplace=True)
    stata2.to_pickle('wiot_'+str(y)+'.pik')

# country list (aside from ROW, which we will construct by aggregating all other countries)
regions = ['USA','CAN','MEX']

# sector aggregation scheme
# A: agriculture
# R: resource extraction
# T: cars + car parts, other transportation equipment
# M: other mfg
# S: services
sectors={'A':[1,2,3],'R':[4,10],'T':[20,21],'M':range(5,10)+range(11,20)+[22],'S':range(23,57)}
#sectors={'A':[1,2,3],'T':[20],'M':range(4,20)+[21,22],'S':range(23,57)}

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

##################################################################################
# aggregate data

intermediates = 0
consumption = 0
investment = 0
value_added = 0
gross_output = 0

print 'Processing data for year...'
for i in range(len(years)):

    y=years[i]
    print '\t'+str(y)


    w = pd.read_pickle('wiot_'+str(y)+'.pik')
    w['row_sector']=w['row_code'].apply(which_sector)
    w['col_sector']=w['col_code'].apply(which_sector)
    w['row_region']=w['row_country'].apply(which_region)
    w['col_region']=w['col_country'].apply(which_region)
    w['col_use']=w['col_code'].apply(which_use)

    # aggregate intermediate inputs by use region/sector and source region/sector...
    # mask for...
    mask = w.col_use=='M' # ...intermediate use
    mask = np.logical_and(mask,w.col_region.isin(regions+['ROW'])) # ... valid use region
    mask = np.logical_and(mask,w.row_sector>=0)  # ...valid source sector
    mask = np.logical_and(mask,w.row_region.isin(regions+['ROW'])) # ... valid source region
    # group by and aggregate
    g = w[mask].groupby(['col_region','col_sector','row_region','row_sector'])
    m = g['value'].sum().reset_index()
    m.rename(columns={'value':'M'},inplace=True)
    m['year'] = y

    if i==0:
        intermediates = m
    else:
        intermediates = intermediates.append(m)
    
    # aggregate consumption by use region and source region/sector...
    # mask for...
    mask = w.col_use=='C' # ...consumption
    mask = np.logical_and(mask,w.col_region.isin(regions+['ROW'])) # ... valid use region
    mask = np.logical_and(mask,w.row_sector>=0)  # ...valid source sector
    mask = np.logical_and(mask,w.row_region.isin(regions+['ROW'])) # ... valid source region
    # group by and aggregate
    g = w[mask].groupby(['col_region','row_region','row_sector'])
    c = g['value'].sum().reset_index()
    c.rename(columns={'value':'C'},inplace=True)
    c['year'] = y

    if i==0:
        consumption = c
    else:
        consumption = consumption.append(c)

    # aggregate investment by use region and source region/sector...
    # mask for...
    mask = w.col_use=='I' # ...investment
    mask = np.logical_and(mask,w.col_region.isin(regions+['ROW'])) # ... valid use region
    mask = np.logical_and(mask,w.row_sector>=0)  # ...valid source sector
    mask = np.logical_and(mask,w.row_region.isin(regions+['ROW'])) # ... valid source region
    # group by and aggregate
    g = w[mask].groupby(['col_region','row_region','row_sector'])
    x = g['value'].sum().reset_index()
    x.rename(columns={'value':'I'},inplace=True)
    x['year'] = y

    if i==0:
        investment = x
    else:
        investment = investment.append(x)
    
    # aggregate VA by use region/sector
    # mask for...
    mask = w.col_use=='M'
    mask = np.logical_and(mask,w.col_region.isin(regions+['ROW'])) # ... valid use region
    mask = np.logical_and(mask,w.row_code.isin([66,67,68,69,70,71]))  # ...valid source sector
    mask = np.logical_and(mask,w.row_region=='TOT') # ... no assigned region
    # group by and aggregate
    g = w[mask].groupby(['col_region','col_sector'])
    va = g['value'].sum().reset_index()
    va.rename(columns={'value':'VA'},inplace=True)
    va['year'] = y

    if i==0:
        value_added = va
    else:
        value_added = value_added.append(va)

    # aggregate GO by use region/sector
    # mask for...
    mask = w.col_use=='M'
    mask = np.logical_and(mask,w.col_region.isin(regions+['ROW'])) # ... valid use region
    mask = np.logical_and(mask,w.row_code==73)  # ...valid source sector
    mask = np.logical_and(mask,w.row_region=='TOT') # ... no assigned region
    # group by and aggregate
    g = w[mask].groupby(['col_region','col_sector'])
    go = g['value'].sum().reset_index()
    go.rename(columns={'value':'GO'},inplace=True)
    go['year'] = y

    if i==0:
        gross_output = go
    else:
        gross_output = gross_output.append(go)

##################################################################################
# merge and check consistency

print 'Checking consistency of aggregated data...'

# merge aggregations with same dimensionality
final_demand = pd.merge(left=consumption,right=investment,
                        how='left',
                        on=['year','col_region','row_region','row_sector'])

output = pd.merge(left=value_added,right=gross_output,
                  how='left',
                  on=['year','col_region','col_sector'])

# market clearing
msums = intermediates.groupby(['year','row_region','row_sector'])['M'].sum().reset_index()
msums.rename(columns={'row_region':'region','row_sector':'sector'},inplace=True)

fsums = final_demand.groupby(['year','row_region','row_sector'])[['C','I']].sum().reset_index()
fsums.rename(columns={'row_region':'region','row_sector':'sector'},inplace=True)

gsums = output[['year','col_region','col_sector','GO']]
gsums = gsums.rename(columns={'col_region':'region','col_sector':'sector'})

sums = pd.merge(left=msums,right=fsums,how='left',on=['year','region','sector'])
sums = pd.merge(left=sums,right=gsums,how='left',on=['year','region','sector'])
sums['diff'] = (sums.GO - sums.M - sums.C - sums.I)
sums['diff'] = sums['diff']/sums['GO']
test = sum(sums['diff']>1e-4)

if test>0:
    print 'Market clearing failure!'

##################################################################################
# save aggregated data

print 'Saving aggregated data to disk...'

intermediates.to_pickle('output/wiod_m.pik')
final_demand.to_pickle('output/wiod_f.pik')
output.to_pickle('output/wiod_vg.pik')
