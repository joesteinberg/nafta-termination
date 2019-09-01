##########################################################################################
# This script does the following:
# - Constructs the IO matrices from processed WIOD data files
# - Writes the matrices to csv files and latex tables
##########################################################################################

##########################################################################################
# Imports, constants, etc.
##########################################################################################
import pandas as pd
import numpy as np
import itertools
import locale
locale.setlocale(locale.LC_ALL,'en_US.utf8')

year=2014
ns=5
sectors=['A','R','T','M','S']
sector_num={'A':0,'R':1,'T':2,'M':3,'S':4}

nc=0
regions=[]
suff = ''

nc=4
regions = {'USA':0,'CAN':1,'MEX':2,'ROW':3}
countries=['USA','CAN','MEX','ROW']

inpath = 'output/'
outpath = 'output/'

m = pd.read_pickle('output/wiod_m.pik')
f = pd.read_pickle('output/wiod_f.pik')
vg = pd.read_pickle('output/wiod_vg.pik')

def assign_region_num(rstr):
    return regions[rstr]

def assign_sector_num(sstr):
    return sector_num[sstr]

##########################################################################################
# load the data
##########################################################################################

vg = vg[vg.year==year]
f = f[f.year==year]
m = m[m.year==year]

vg['col_region_num'] = vg['col_region'].apply(assign_region_num)
vg['col_sector_num'] = vg['col_sector'].apply(assign_sector_num)

f['col_region_num'] = f['col_region'].apply(assign_region_num)
f['row_region_num'] = f['row_region'].apply(assign_region_num)
f['row_sector_num'] = f['row_sector'].apply(assign_sector_num)

m['col_region_num'] = m['col_region'].apply(assign_region_num)
m['row_region_num'] = m['row_region'].apply(assign_region_num)
m['col_sector_num'] = m['col_sector'].apply(assign_sector_num)
m['row_sector_num'] = m['row_sector'].apply(assign_sector_num)

# set investment demand for agriculture and consumption of resources to zero
f.loc[f.row_sector=='A','I'] = 0.0
#f.loc[f.row_sector=='R','C'] = 0.0

##########################################################################################
# construct the input-output matrix
##########################################################################################

vg = vg[vg.year==year].sort_values(by=['col_region_num','col_sector_num']).reset_index()
f = f[f.year==year].sort_values(by=['col_region_num','row_region_num','row_sector_num']).reset_index()
m = m[m.year==year].sort_values(by=['col_region_num','col_sector_num','row_region_num','row_sector_num']).reset_index()

rowsums = np.zeros( nc*ns + 1 )
colsums = np.zeros( nc*ns + nc*2 )

MM = m.pivot_table(values='M', index=['row_region_num','row_sector_num'], columns=['col_region_num','col_sector_num'])
VV = vg['VA'].values.reshape((1,nc*ns))
FF = f.pivot_table(values=['C','I'], index=['row_region_num','row_sector_num'], columns=['col_region_num'])
VV = np.hstack((VV,np.zeros((1,nc*2))))

iomat=np.vstack( ( np.hstack((MM,FF)) , VV ) )

for row in range(0,nc*ns + 1):
    rowsums[row] = np.sum(iomat[row,:])

for col in range(0,nc*ns + nc*2):
    colsums[col] = np.sum(iomat[:,col])

##########################################################################################
# rebalance
##########################################################################################

def coeffs(iomat):
    # Given world IO matrix (iomat), calculates IO coefficients and returs them in A

    A=np.zeros(iomat.shape)
    for col in range(0,A.shape[1]):
        A[:,col] = iomat[:,col]/np.sum(iomat[:,col])
    return A

def ras(iomat0,rowsums1,colsums1):
    # Given an initial IO matrix (iomat), and desired rowsums (rowsums1) and colsums (colsums1),
    # performs the RAS balancing procedure. Returns a new IO matrix (iomat) that is consistent
    # with desired row- and colsums.

    A0 = coeffs(iomat0)
    iomat = np.dot(A0,np.diag(colsums1))

    go=True
    iter=0
    maxit=10000
    tol=1.0e-8

    while go:
        iter=iter+1
        rowsums = np.sum(iomat,axis=1)
        r = np.divide(rowsums1,rowsums)
        iomat = np.dot(np.diag(r),iomat)
        colsums = np.sum(iomat,axis=0)
        s = np.divide(colsums1,colsums)
        iomat = np.dot(iomat,np.diag(s))
        colsums = np.sum(iomat,axis=0)
        rowsums = np.sum(iomat,axis=1)

        norm1 = max(np.divide(abs(rowsums-rowsums1),rowsums1))
        norm2 = max(np.divide(abs(colsums-colsums1),colsums1))
        if((norm1 <tol and norm2 <tol) or iter == maxit):
            go=False

    if iter==maxit:
        print 'RAS iteration did not converge!'
        print 'iter = ', iter, ' diff = ', max(norm1,norm2)
    else:
        print 'RAS converged after ',str(iter),' iterations'


    return iomat

# make sure it's balanced after imposing the restrictions on construction usage...
colsums[0:(nc*ns)] = rowsums[0:(nc*ns)] # make sure markets clear: gross output = total demand for each country/sector
rowsums[-1] = colsums[(nc*ns):].sum() # world value added must equal world final demand
iomat2 = ras(iomat,rowsums,colsums) # run RAS

##########################################################################################
# write output
##########################################################################################

def write_iomat_csv(iomat,fname):
    # Write world IO matrix (iomat) to csv file called filename

    usgdp = iomat[-1,0:ns].sum()
    iomat2 = np.vstack((iomat,np.sum(iomat,axis=0).reshape((1,nc*ns+nc*2))))
    iomat2 = np.hstack((iomat2,np.sum(iomat2,axis=1).reshape((nc*ns+2,1))))
    iomat2 = 100*iomat2/usgdp
    np.savetxt(fname=outpath+fname,X=iomat2,fmt='%0.15f',delimiter=' ')

def write_iomat_latex(iomat,rowsums,colsums,caption,units,label,fname):
    # Given a world IO matrix (iomat), rowsums, colsums, creates a latex file
    # in location filename that contains a table with given caption and label.

    usgdp = iomat[-1,0:ns].sum()
    iomat2 = 100*iomat[:,:]/usgdp
    rowsums2 = 100*rowsums/usgdp
    colsums2 = 100*colsums/usgdp

    M=iomat2[0:(nc*ns),0:(nc*ns)]
    V=iomat2[-1,0:(nc*ns)]
    Fc=iomat2[0:(nc*ns),(nc*ns):+((nc*ns)+nc)]
    Fx=iomat2[0:(nc*ns),((nc*ns)+nc):]

    with open(outpath + fname + '_m.tex','wb') as file:
        file.write('\\begin{landscape}\n')
        file.write('\\begin{table}[p]\n')
        #file.write('\\renewcommand{\\arraystretch}{1.2}\n')
        file.write('\\begin{center}\n')
        file.write('\\caption{'+caption+', intermediate inputs portion '+units+'}\n')
        file.write('\\label{tab:'+label+'_m}\n')

        file.write('\\footnotesize\n')
        file.write('\\begin{tabular}{cc')
        for i in range(0,nc*ns):
            file.write('c')
        file.write('}\n')
        file.write('\\toprule\n')
            
        # country names
        file.write('&')
        for c in countries:
            file.write('& \\multicolumn{'+str(ns)+'}{c}{'+c+'}')
        #for c in countries:
        #    file.write('& \\multicolumn{2}{c}{'+c+'}')
        file.write('\\\\\n')
        
        # underline country names
        x=3
        for i in range(nc):
            file.write('\\cmidrule(rl){'+str(x)+'-'+str(x+ns-1)+'}')
            x=x+ns
        #for i in range(nc):
        #    file.write('\\cmidrule(rl){'+str(x)+'-'+str(x+1)+'}')
        #    x=x+2
        file.write('\n')

        # sector names
        file.write('&')
        for c in countries:
            for s in sectors:
                file.write('&' + s)
        #for c in countries:
        #    file.write('& Cons & Inv')
        file.write('\\\\\n')
        file.write('\\midrule\n')

        for i in range(0,nc):
            file.write('\\multirow{'+str(ns)+'}{*}{\\begin{sideways}'+countries[i]+'\\end{sideways}}')
            for ii in range(0,ns):
                file.write('&'+sectors[ii])
                for j in range(0,nc):
                    for jj in range(0,ns):
                        tmpstr = '-'
                        if M[i*ns+ii][j*ns+jj] > 1e-6:
                            tmpstr = locale.format('%0.2f',M[i*ns+ii,j*ns+jj],grouping=True)
                        file.write('&'+tmpstr)
                #for j in range(0,nc):
                #    tmpstr = locale.format('%0.2f',Fc[i*ns+ii,j],grouping=True)
                #    file.write('&'+tmpstr)
                #    tmpstr = locale.format('%0.2f',Fx[i*ns+ii,j],grouping=True)
                #    file.write('&'+tmpstr)
                file.write('\\\\\n')
            file.write('\\midrule\n')
            
        #file.write('\\midrule\n')
        file.write('VA &')
        for i in range(0,nc):
            for ii in range(0,ns):
                tmpstr='-'
                if V[i*ns+ii]>1e-6:
                    tmpstr = locale.format('%0.2f',V[i*ns+ii],grouping=True)
                file.write('&'+tmpstr)

        file.write('\\\\\n')

        file.write('\\midrule\n')
        file.write('GO &')
        for i in range(0,nc):
            for ii in range(0,ns):
                tmpstr='-'
                if colsums2[i*ns+ii]>1e-6:
                    tmpstr = locale.format('%0.2f',colsums2[i*ns+ii],grouping=True)
                file.write('&'+tmpstr)

        file.write('\\\\\n')
            
        file.write('\\bottomrule\n')
        file.write('\\end{tabular}\n')
        file.write('\\end{center}\n')
        file.write('\\end{table}\n')
        file.write('\\end{landscape}\n')

    with open(outpath + fname + '_fc.tex','wb') as file:
        #file.write('\\begin{landscape}\n')
        file.write('\\begin{table}[p]\n')
        #file.write('\\renewcommand{\\arraystretch}{1.2}\n')
        file.write('\\begin{center}\n')
        file.write('\\caption{'+caption+', final demand portion '+units+'}\n')
        file.write('\\label{tab:'+label+'_fc}\n')

        file.write('\\footnotesize\n')
        file.write('\\begin{tabular}{cc')
        for i in range(0,nc*2):
            file.write('c')
        file.write('}\n')
        file.write('\\toprule\n')
            
        # country names
        file.write('&')
        for c in countries:
            file.write('& \\multicolumn{2}{c}{'+c+'}')
        file.write('\\\\\n')
        
        # underline country names
        x=3
        for i in range(nc):
            file.write('\\cmidrule(rl){'+str(x)+'-'+str(x+1)+'}')
            x=x+2
        file.write('\n')

        # sector names
        file.write('&')
        for c in countries:
            file.write('& Cons & Inv')
        file.write('\\\\\n')
        file.write('\\midrule\n')

        for i in range(0,nc):
            file.write('\\multirow{'+str(ns)+'}{*}{\\begin{sideways}'+countries[i]+'\\end{sideways}}')
            for ii in range(0,ns):
                file.write('&'+sectors[ii])
                for j in range(0,nc):

                    tmpstr='-'
                    if Fc[i*ns+ii,j]>1e-6:
                        tmpstr = locale.format('%0.2f',Fc[i*ns+ii,j],grouping=True)
                    file.write('&'+tmpstr)

                    tmpstr='-'
                    if Fx[i*ns+ii,j]>1e-6:
                        tmpstr = locale.format('%0.2f',Fx[i*ns+ii,j],grouping=True)
                    file.write('&'+tmpstr)
                file.write('\\\\\n')
            file.write('\\midrule\n')
            
        #file.write('\\midrule\n')
        file.write('VA &')
        for i in range(0,nc):
            file.write('& - & -')
        file.write('\\\\\n')

        file.write('\\midrule\n')
        file.write('GO &')
        for i in range(0,nc):
            file.write('& - & -')

        file.write('\\\\\n')
            
        file.write('\\bottomrule\n')
        file.write('\\end{tabular}\n')
        file.write('\\end{center}\n')
        file.write('\\end{table}\n')
        #file.write('\\end{landscape}\n')
        file.write('\\normalsize\n')

if(year==2000):
    write_iomat_csv(iomat2,'iomat_old.txt')
else:
    write_iomat_csv(iomat2,'iomat.txt')
    write_iomat_latex(iomat2,
                      rowsums,
                      colsums,
                      '2014 input-output table',
                      '(U.S. GDP = 100)',
                      'iomat',
                      'iomat')


##########################################################################################
# alternative no-IO matrix
##########################################################################################

if year!=2000:
    m.loc[:,'M']=0.0

    vg = vg[vg.year==year].sort_values(by=['col_region_num','col_sector']).reset_index()
    f = f[f.year==year].sort_values(by=['col_region_num','row_region_num','row_sector']).reset_index()
    m = m[m.year==year].sort_values(by=['col_region_num','col_sector','row_region_num','row_sector']).reset_index()

    vg['col_region_num'] = vg['col_region'].apply(assign_region_num)
    f['col_region_num'] = f['col_region'].apply(assign_region_num)
    f['row_region_num'] = f['row_region'].apply(assign_region_num)
    m['col_region_num'] = m['col_region'].apply(assign_region_num)
    m['row_region_num'] = m['row_region'].apply(assign_region_num)
    
    rowsums = np.zeros( nc*ns + 1 )
    colsums = np.zeros( nc*ns + nc*2 )

    MM = m.pivot_table(values='M', index=['row_region_num','row_sector'], columns=['col_region_num','col_sector'])
    VV = vg['VA'].values.reshape((1,nc*ns))
    FF = f.pivot_table(values=['C','I'], index=['row_region_num','row_sector'], columns=['col_region_num'])
    VV = np.hstack((VV,np.zeros((1,nc*2))))

    iomat=np.vstack( ( np.hstack((MM,FF)) , VV ) )

    for row in range(0,nc*ns + 1):
        rowsums[row] = np.sum(iomat[row,:])

    for col in range(0,nc*ns + nc*2):
        colsums[col] = np.sum(iomat[:,col])

    # make sure it's balanced after imposing the restrictions on construction usage...
    colsums[0:(nc*ns)] = rowsums[0:(nc*ns)] # make sure markets clear
    rowsums[-1] = colsums[(nc*ns):].sum() # world value added must equal world final demand
    iomat2 = ras(iomat,rowsums,colsums) # run RAS

    write_iomat_csv(iomat2,'iomat_noio.txt')
