#!/usr/bin/env python
# coding: utf-8

# In[102]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import glob
import re


# In[2]:


def read_files(path):
    try:
        freqs = np.array(open(path+'freqs.txt','rb').read().split(',')).astype(float)
    except:
        print "no freqs file found"
        return ""
    data_files = glob.glob(path+'data*.txt')
    if len(data_files)== 0:
        print "no data files"
        return ""
    data = [0]*(len(data_files))
    for f in data_files:
        ind = int(re.findall('.*data.?_(\d+).txt',f)[0])
        temp = np.array(open(f,'rb').read().split(',')).astype(float)
        temp = [temp[i]+1j*temp[i+1] for i in xrange(0,len(temp),2)]
        data[ind] = np.array(temp) 
    res = [freqs] + data
    return res


# In[64]:


path = 'N3250A/meas/12800 10cell/'
data_files = glob.glob(path+'data*.txt')
r = read_files(path)


# In[66]:


np.abs(r[1:])


# In[70]:


freqs = r[0]
vols = r[1:]
ls = [1.53e-3]*33
ls =np.round(np.cumsum(ls),5)


# In[104]:


plt.close('all')
plt.figure(figsize = (10,10))
X,Y = np.meshgrid(freqs,ls)
plt.contourf(X,Y,np.abs(vols),locator=ticker.LogLocator(base = 2),cmap = plt.cm.YlGnBu_r)
plt.xlabel("frequency")
plt.ylabel('distance')
plt.colorbar()
plt.show()


# In[98]:


plt.close('all')
#plt.figure(figsize=[10,10])
#plt.plot(vs_mat[1,:])
#plt.show()


# In[90]:


vs_mat = np.matrix(np.abs(vols))


# In[97]:


vs_mat[1,:]


# In[86]:


get_ipython().run_line_magic('matplotlib', 'notebook')
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure(figsize=[10,10])

ax = fig.gca(projection='3d')
vols = np.array(vols)
ax.plot_surface(X, Y, np.abs(vols), cmap = plt.cm.YlGnBu_r)


plt.show()

