import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.mlab as mlab
from matplotlib import colors
from matplotlib.pyplot import MultipleLocator


name1='BTE.kappa_glass_resolved'
data1= pd.read_csv(name1,sep='\s+',header=None)

data1.columns=['w1','w2','v1','v2','v3','v4','v5','v6','v7','v8','v9']
data1['v9_sum'] = data1.groupby(['w1','w2'])['v9'].transform('sum')
#data1['avg'] = (data1['v1']+ data1['v5']+ data1['v9'])/3
data1['avg'] = data1['v1']
data2 = data1[data1['avg']>10**(-10)]
# Frequency unit: THz
plt.scatter(data2['w1']*0.15915494327376,data2['w2']*0.15915494327376,c=data2['avg'],
norm=colors.LogNorm(vmax=10**(-4),vmin=10**(-10)),cmap='jet',s=0.1)
plt.colorbar()
ax1=plt.gca()
ax1.spines['bottom'].set_linewidth(1)
ax1.spines['top'].set_linewidth(1)
ax1.spines['right'].set_linewidth(1)
ax1.spines['left'].set_linewidth(1)

major = MultipleLocator(2)  
minor = MultipleLocator(1) 


ax1.xaxis.set_minor_locator(minor)
ax1.xaxis.set_major_locator(major)

ax1.yaxis.set_minor_locator(minor)
ax1.yaxis.set_major_locator(major)
ax1.tick_params(direction='in', length=2, width=1, colors='k',which='minor')
ax1.tick_params(direction='in', length=4, width=1, colors='k',which='major')


#plt.xlim(-0.5,13.5)
#plt.ylim(-0.5,13.5)
#plt.xticks=([])
#plt.yticks=([])
#ax1.set_xticks([])
#ax1.set_yticks([])
plt.show()
#data1['v1_sum'] = data1.groupby(['w1','w2'])['v1'].transform('sum')
#df0 = data1.groupby(['w1','w2']).agg("sum")
#print(data1.head())
