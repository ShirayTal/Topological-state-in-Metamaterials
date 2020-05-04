#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from winsound import Beep


# In[2]:


get_ipython().run_line_magic('run', 'parse_cti.ipynb')
get_ipython().run_line_magic('run', 'metamaterials_funs.ipynb')


# In[3]:


#w=100
R1 = 49.0357
C1 = 3.51627e-11
L1 = 7.32911e-07
G1 = 0.005312


#w=750
R2 = 3.86443
C2 = 8.3313e-11
L2 = 3.16413e-07
G2 = 0.0056311

#sma cable
C3 = 0.96E-12
L3 = 50**2*0.96E-12 

FREQS = np.arange(0.1,20,0.001)*10**9


# In[64]:


def calc_z0(omega,r,g,c,l):
    return np.sqrt((r+1j*omega*l)/(g+1j*omega*c))
    #return np.sqrt(l/c)

def calc_prop_coeff(omega,r,g,c,l):
    return np.sqrt((r+1j*omega*l)*(g+1j*omega*c))
    #return 1j*omega*np.sqrt(l*c)

def calc_abcd(z0, beta=0, dist=0, comp ="TL", omega = 1):
    # calc abcd matrix from z0, d and beta
    if comp == "Sr":
        a = 1
        b = z0
        c = 0
        d = 1
    if comp == "Sh":
        a = 1
        b = 0
        c = 1/z0
        d = 1
    if comp == "R":
        a = 1
        b = z0
        c = 0
        d = 1
    if comp == "C":
        a = 1
        b = 1/(1j*omega*z0)
        c = 0
        d = 1
    if comp == "L":
        a = 1
        b = 0
        c = 1j*omega*z0
        d = 1
    if comp == "TL":
        a = np.cosh(beta*dist)
        b = z0*np.sinh(beta*dist)
        c = np.sinh(beta*dist)/z0
        d = np.cosh(beta*dist)
    abcd = np.matrix([[a,b],[c,d]])
    return abcd

def abcd_to_z(abcd):
    # calc Z matrix from abcd matrix
    z = np.matrix([[1j,1j],[1j,1j]])
    c = abcd[1,0]
    z[0,0] = abcd[0,0]/c
    z[0,1] = np.linalg.det(abcd)/c
    z[1,0] = 1/c
    z[1,1] = abcd[1,1]/c
    return z

def z_to_s(z):
    # calc S matrix from Z matrix
    u = np.eye(2,2)
    temp = z + u
    s = np.matmul(temp.getI(),z-u)
    return s

def z_to_abcd(z):
    # calc abcd matrix from z
    abcd = np.matrix([[1j,1j],[1j,1j]])
    z21 = z[1,0]
    abcd[0,0]=z[0,0]/z21
    abcd[0,1] = np.linalg.det(z)/z21
    abcd[1,0] = 1/z21
    abcd[1,1] = abcd[1,1]/z21
    return abcd

def CPW_total_abcd(types = [], sequence = []):
    tot_abcd = np.eye(2,2)
    for i in sequence:
        cur_abcd = types[i]
        tot_abcd = np.dot(tot_abcd,cur_abcd)
    return tot_abcd

def generate_abcd_chain(omega, gammas, zs,ls, comps):
    abcd_chain = [np.eye(2,2)]
    for i in xrange(len(gammas)):
        cur_mat = calc_abcd(zs[i],gammas[i],ls[i],comps[i],omega)
        abcd_chain.append(np.dot(abcd_chain[i],cur_mat))
    return abcd_chain[1:]

def get_vr(abcd, z_in,z_l):
    vec = np.dot(abcd, [1,1/z_l])
    v = vec[0,0]
    i = vec[0,1]
    v_in = 0.5*(v+z_in*i)
    vr = z_l/v_in
    return vr

def generate_Vchain(omega, gammas, zs, ls, comps, v_in):
    z_out = 50    
    abcd_chain = generate_abcd_chain(omega, gammas, zs, ls, comps)
    v_r = get_vr(abcd_chain[-1],50,z_out)
    v_out = v_in*v_r
    v_chain = [v_out]
    out_vec = [v_out, v_out/z_out]
    for mat in abcd_chain:
        cur_vec = np.dot(mat,out_vec)
        v_chain.append(cur_vec[0,0])
    return v_chain[-1::-1]

def calc_coeff(abcd):
    [v_start,i_start] = get_start_params(abcd)
    t = 2/(v_start+50*i_start)
    r = (v_start-50*i_start)/(v_start+50*i_start)
    return r,t

def get_start_params(abcd):
    e = np.array((1,1/50.0))
    [v_start,i_start] = np.abs(np.array(np.dot(abcd,e))[0])
    return [v_start,i_start]

def generate_zChain(types, dists, ls):
    edge = 1
    z_chain =  []
    for dist in dists:
        pass


# # T and R of 10 cell

# In[7]:


freqs = FREQS
d = 2.6*10**-3
rs = []
ts = []
for freq in freqs:
    omega = 2*np.pi*freq
    Z1 = calc_z0(omega, R1, G1, C1, L1)
    beta1 = calc_prop_coeff(omega, R1, G1, C1, L1)
    nar1 = calc_abcd(Z1, beta1, d)
    Z2 = calc_z0(omega, R2, G2, C2, L2)
    beta2 = calc_prop_coeff(omega, R2, G2, C2, L2)
    wide1 = calc_abcd(Z2, beta2, d)

    types = [nar1, wide1]
    ten_cell = [0,1]*10
    a = CPW_total_abcd(types, ten_cell)
    cur_r,cur_t = calc_coeff(a)
    rs.append(cur_r)
    ts.append(cur_t)


# In[8]:


x,y = parse_cti('10cell_first.cti')
y = complexify_y(y)
r = y[0]
t = y[1]


# In[9]:


plt.close('all')
plt.figure(figsize = [10,10])
plt.scatter(freqs,np.abs(ts), label = 'simulation')
plt.scatter(x, np.abs(t), label = 'experiment')
plt.title('10 cell transfer')
plt.legend()
plt.xlabel('frequency [Hz]')
plt.ylabel('|t|')
plt.show()


# In[56]:


plt.close('all')
plt.figure(figsize = [10,10])
plt.scatter(freqs,np.abs(rs), label = 'simulation')
plt.scatter(x, np.abs(r), label = 'experiment' )
plt.title('10 cell return')
plt.legend()
plt.xlabel('frequency [Hz]')
plt.ylabel('|r|')
plt.show()


# ## T and R of 4 cell

# In[57]:


freqs = FREQS
d = 6.6*10**-3
rs = []
ts = []
for freq in freqs:
    omega = 2*np.pi*freq
    Z1 = calc_z0(omega, R1, G1, C1, L1)
    beta1 = calc_prop_coeff(omega, R1, G1, C1, L1)
    nar1 = calc_abcd(Z1, beta1, d)
    Z2 = calc_z0(omega, R2, G2, C2, L2)
    beta2 = calc_prop_coeff(omega, R2, G2, C2, L2)
    wide1 = calc_abcd(Z2, beta2, d)

    types = [nar1, wide1]
    four_cell = [0,1]*4
    a = CPW_total_abcd(types, four_cell)
    cur_r,cur_t = calc_coeff(a)
    rs.append(cur_r)
    ts.append(cur_t)


# In[58]:


x,y = parse_cti('4cell_first.cti')
y = complexify_y(y)
r = y[0]
t = y[1]


# In[59]:


plt.close('all')
plt.figure(figsize = [10,10])
plt.scatter(freqs,np.abs(ts), label = 'simulation')
plt.scatter(x, np.abs(t), label = 'experiment')
plt.title('4 cell transfer')
plt.legend()
plt.xlabel('frequency [Hz]')
plt.ylabel('|t|')
plt.show()


# In[60]:


plt.close('all')
plt.figure(figsize = [10,10])
plt.scatter(freqs,np.abs(rs), label = 'simulation')
plt.scatter(x, np.abs(r), label = 'experiment' )
plt.title('4 cell return')
plt.legend()
plt.xlabel('frequency [Hz]')
plt.ylabel('|r|')
plt.show()


# # T and R of 4-10 cell

# In[15]:


freqs = FREQS
d1 = 6.6*10**-3
d2 = 2.6*10**-3
rs = []
ts = []
for freq in freqs:
    omega = 2*np.pi*freq
    Z1 = calc_z0(omega, R1, G1, C1, L1)
    beta1 = calc_prop_coeff(omega, R1, G1, C1, L1)
    nar1 = calc_abcd(Z1, beta1, d1)
    nar2 = calc_abcd(Z1, beta1, d2)
    Z2 = calc_z0(omega, R2, G2, C2, L2)
    beta2 = calc_prop_coeff(omega, R2, G2, C2, L2)
    wide1 = calc_abcd(Z2, beta2, d1)
    wide2 = calc_abcd(Z2, beta2, d2)
    Z3 = 50
    beta3 = calc_prop_coeff(omega,0,0,C3,L3)
    d3 = 76e-3
    wire = calc_abcd(Z3, beta3,d3)
    types = [nar1, wide1, nar2,wide2,wire]
    four_ten_cell = [0,1]*4+[4]+[2,3]*10
    a = CPW_total_abcd(types, four_ten_cell)
    cur_r,cur_t = calc_coeff(a)
    rs.append(cur_r)
    ts.append(cur_t)


# In[11]:


x,y = parse_cti('4_10_cell.cti')
y = complexify_y(y)
r = y[0]
t = y[1]


# In[14]:


plt.close('all')
plt.figure(figsize = [10,10])
plt.scatter(freqs,np.abs(ts), label = 'simulation')
plt.scatter(x, np.abs(t), label = 'experiment')
plt.title('4-10 cell transfer')
plt.legend()
plt.xlabel('frequency [Hz]')
plt.ylabel('|t|')
plt.show()


# In[16]:


plt.close('all')
plt.figure(figsize = [10,10])
plt.scatter(freqs,np.abs(rs), label = 'simulation')
plt.scatter(x, np.abs(r), label = 'experiment' )
plt.title('4 10 cell return')
plt.legend()
plt.xlabel('frequency [Hz]')
plt.ylabel('|r|')
plt.show()


# # T and R of 10-4 cell

# In[17]:


freqs = FREQS
d1 = 6.6*10**-3
d2 = 2.6*10**-3
rs = []
ts = []
for freq in freqs:
    omega = 2*np.pi*freq
    Z1 = calc_z0(omega, R1, G1, C1, L1)
    beta1 = calc_prop_coeff(omega, R1, G1, C1, L1)
    nar1 = calc_abcd(Z1, beta1, d1)
    nar2 = calc_abcd(Z1, beta1, d2)
    Z2 = calc_z0(omega, R2, G2, C2, L2)
    beta2 = calc_prop_coeff(omega, R2, G2, C2, L2)
    wide1 = calc_abcd(Z2, beta2, d1)
    wide2 = calc_abcd(Z2, beta2, d2)
    Z3 = 50
    beta3 = calc_prop_coeff(omega,0,0,C3,L3)
    d3 = 76e-3
    wire = calc_abcd(Z3, beta3,d3)
    types = [nar1, wide1, nar2,wide2,wire]
    four_ten_cell = [2,3]*10+[4]+[0,1]*4
    a = CPW_total_abcd(types, four_ten_cell)
    cur_r,cur_t = calc_coeff(a)
    rs.append(cur_r)
    ts.append(cur_t)


# In[18]:


x,y = parse_cti('10_4_cell.cti')
y = complexify_y(y)


# In[19]:


plt.close('all')
plt.figure(figsize = [10,10])
plt.scatter(freqs,np.abs(ts), label = 'simulation')
plt.scatter(x, np.abs(t), label = 'experiment')
plt.title('10-4 cell transfer')
plt.legend()
plt.xlabel('frequency [Hz]')
plt.ylabel('|t|')
plt.show()


# In[20]:


plt.close('all')
plt.figure(figsize = [10,10])
plt.scatter(freqs,np.abs(rs), label = 'simulation')
plt.scatter(x, np.abs(r), label = 'experiment' )
plt.title('10 4 cell return')
plt.legend()
plt.xlabel('frequency [Hz]')
plt.ylabel('|r|')
plt.show()


# # Sweep of 10 cell

# In[129]:


vs = []
omegas = np.arange(0.1,2,0.00015)*1e10
for omega in omegas:
    omega = 2*np.pi*omega
    Z1 = calc_z0(omega, R1, G1, C1, L1)
    beta1 = calc_prop_coeff(omega, R1, G1, C1, L1)
    Z2 = calc_z0(omega, R2, G2, C2, L2)
    beta2 = calc_prop_coeff(omega, R2, G2, C2, L2)
    gammas = [beta1,beta1,beta1,beta1,beta2,beta2,beta2,beta2]*10
    ls = [(1.3/2)*1e-3]*80
    zs = [Z1,Z1,Z1,Z1,Z2,Z2,Z2,Z2]*10
    comp = ["TL"]*80
    v = np.abs(generate_Vchain(omega,gammas,zs,ls,comp,1)[1:])
    vs.append(v)


# In[130]:


len(vs[0])


# In[131]:


plt.close('all')
ls = [(1.3/2)*1e-3]*80
plt.figure(figsize = (15,15))
Y,X= np.meshgrid(np.cumsum(ls)*1e3, omegas)
plt.contourf(X,Y,vs,locator=ticker.LogLocator(base = 2), cmap = plt.cm.RdBu)

plt.colorbar()
plt.show()


# In[35]:


plt.close('all')
plt.plot(np.cumsum(ls)*1e3,vs[1000])
plt.show()


# In[54]:


plt.close('all')
plt.scatter(range(len(vs[1000])),vs[1000])
plt.show()


# In[82]:


omega = 0.25e10
Z1 = calc_z0(omega, R1, G1, C1, L1)
beta1 = calc_prop_coeff(omega, R1, G1, C1, L1)
Z2 = calc_z0(omega, R2, G2, C2, L2)
beta2 = calc_prop_coeff(omega, R2, G2, C2, L2)
gammas = [beta1,beta2]*4
ls = [d1,d2]*4
zs = [Z1,Z2]*4
comp = ["TL"]*8
v = np.abs(generate_Vchain(omega,gammas,zs,ls,comp,1))


# # Testing defects

# In[229]:


#new cpw
# w = 50 s = 1000
R = 94.4763
L = 7.75781e-07
C = 3.25918e-11
G = 0.00494502


# In[307]:


get_ipython().run_line_magic('time', '')
freqs = np.arange(0.1,2,0.00025)*10**10
d = 10*10**-3
rs = {'ref':[],'def':[]}
ts = {'ref':[],'def':[]}
for freq in freqs:
    omega = 2*np.pi*freq
    Z1 = calc_z0(omega, R1, G1, C1, L1)
    beta1 = calc_prop_coeff(omega, R1, G1, C1, L1)
    nar1 = calc_abcd(Z1, beta1, d)
    Z2 = calc_z0(omega, R2, G2, C2, L2)
    beta2 = calc_prop_coeff(omega, R2, G2, C2, L2)
    wide1 = calc_abcd(Z2, beta2, d)
    Z3 = calc_z0(omega,R, G, C, L)
    beta3 = calc_prop_coeff(omega,R, G, C, L)
    add = calc_abcd(Z3,beta3,d,comp = 'L', omega=omega)
    #add = np.eye(2,2)
    types = [nar1, wide1, add]
    #ref = [1,0]*2+[0,1]*6+[1,0]*2
    ref = [1,0]*10
    a = CPW_total_abcd(types, ref)
    ref_r,ref_t = calc_coeff(a)
    rs['ref'].append(ref_r)
    ts['ref'].append(ref_t)
    
    #defect = [1,0]*2+[0,1]*6+[1,0]*2
    #defect = [1,0]*4+[1,0,1,2]+[1,0]*4
    #defect = [0,1,2]+[0,1,2]+[0,1,2]+[0,1,2]
    defect = [1,0]*9+[1,2]
    b = CPW_total_abcd(types, defect)
    def_r,def_t = calc_coeff(b)
    rs['def'].append(def_r)
    ts['def'].append(def_t)
Beep(880,200)
Beep(660,200)


# In[308]:


plt.close('all')
plt.figure(figsize = [20,10])
plt.scatter(freqs,np.abs(ts['ref']), label = 'ref')
plt.scatter(freqs,np.abs(ts['def']), label = 'defect')
#plt.scatter(x, np.abs(t), label = 'experiment')
plt.title('reference vs defect')
plt.legend(fontsize = 20)
plt.xlabel('frequency [Hz]')
plt.ylabel('|t|')
plt.xticks(np.arange(0,2,0.1)*10**10)
plt.show()


# In[85]:


plt.close('all')
plt.figure(figsize = [20,10])
plt.scatter(freqs,np.abs(rs['ref']), label = 'ref')
plt.scatter(freqs,np.abs(rs['def']), label = 'defect')
#plt.scatter(x, np.abs(t), label = 'experiment')
plt.title('12 8mm')
plt.legend()
plt.xlabel('frequency [Hz]')
plt.ylabel('|r|')
plt.show()


# In[265]:


plt.close('all')
plt.figure(figsize = [20,10])
plt.scatter(freqs,np.abs(ts), label = 'simulation')
#plt.scatter(x, np.abs(t), label = 'experiment')
plt.title('10 1cm')
plt.legend()
plt.xlabel('frequency [Hz]')
plt.ylabel('|t|')
plt.show()


# # One mass Metamaterial with defect

# In[179]:


# one mass metamaterial
Ls = 1.89276e-07
Rs = 10
Cs = 2.02895e-10
Gs = 0.0308397
L = 15e-9
C = 1e-12
C_n = 2.5e-9 # the defect (need to chane the name)
d = 11e-3

rs2 = []
ts2 = []
freqs = np.arange(0.1,10,0.001)*10**9
for freq in freqs:
    omega = 2*np.pi*freq
    Zs = calc_z0(omega,Rs,Gs,Cs,Ls)
    beta_s = calc_prop_coeff(omega,Rs,Gs,Cs,Ls)
    TL = calc_abcd(Zs,beta_s,d/2.0,'TL',omega)
    Coil = calc_abcd(1j*omega*L,comp = 'Sh',omega=omega)
    Cap = calc_abcd(1/(1j*omega*C),comp = 'Sr',omega=omega)
    Def = calc_abcd(1j*omega*C_n,comp = 'Sh',omega=omega)
    types = [TL, Coil, Cap, Def]
    CLRH = [2,0,1,1,0]*4+[2,0,3,3,0]*1+[2,0,1,1,0]*3
    a = CPW_total_abcd(types,CLRH)
    r,t = calc_coeff(a)
    rs2.append(r)
    ts2.append(t)


# In[287]:


vs_clean = []
omegas = np.arange(0.001,1,0.00015)*1e10
for omega in omegas:
    omega = 2*np.pi*omega
    Zc = 1/(1j*omega*C)
    D = 1j*omega*C_n
    Zl = 1j*omega*L
    Z0 = calc_z0(omega,Rs,Gs,Cs,Ls)
    gamma = calc_prop_coeff(omega,Rs,Gs,Cs,Ls)
    gammas = [0,gamma,0,0,gamma]*8
    ls = [0,5.5e-3,0,0,5.5e-3]*8
    zs = [Zc,Z0,Zl,Zl,Z0]*4+[Zc,Z0,D,D,Z0]*0+[Zc,Z0,Zl,Zl,Z0]*4
    comp = ["Sr","TL","Sh","Sh","TL"]*8
    v = np.abs(generate_Vchain(omega,gammas,zs,ls,comp,1)[1:])
    vs_clean.append(choose_cell(v,[0,4],5))
    


# In[297]:


plt.close('all')
plt.figure(figsize = (15,15))
Y,X= np.meshgrid(np.cumsum(choose_cell(ls,[1,4],5))*1e3, omegas)
plt.contourf(X,Y,vs_clean,locator=ticker.LogLocator(base = 20), cmap = plt.cm.RdBu_r)

plt.colorbar()
plt.show()


# In[306]:


plt.close('all')
plt.figure(figsize = (15,15))
Y,X= np.meshgrid(np.cumsum(choose_cell(ls,[1,4],5))*1e3, omegas)
plt.contourf(X,Y,vs,locator=ticker.LogLocator(base = 20), cmap = plt.cm.RdBu)

plt.colorbar()
plt.show()


# In[299]:



for i in xrange(len(choose_cell(ls,[1,4],5))):
    plt.close('all')
    #plt.plot(np.array(np.matrix(vs)[:,9]), label = 'defect')
    plt.plot(np.array(np.matrix(vs)[:,i]), label = 'other')
    plt.yscale('log')
    plt.legend()
    plt.savefig('def/im'+str(i)+".jpg")


# In[304]:


plt.close('all')
plt.plot(np.array(np.matrix(vs)[1977,:])[0])
plt.yscale('log')
plt.show()


# In[305]:


plt.close('all')
plt.figure(figsize = [10,10])
plt.scatter(freqs,np.abs(ts), label = 'clean')
#plt.scatter(freqs,np.abs(ts1), label = 'def1')

plt.scatter(freqs,np.abs(ts2), label = 'def2')
plt.legend()
plt.show()


# In[155]:


popt, pcov = opt.curve_fit(fourier1,range(80),vs[20],p0 = [100, 100,0,0.01], method ='trf' )
fit = [fourier1(x, *popt) for x in xrange(80)]
plt.close('all')
plt.plot(np.array(np.matrix(vs)[20,:])[0])
plt.plot(fit)
#plt.yscale('log')
plt.show()


# In[153]:


plt.plot(np.array(np.matrix(vs)[20,:])[0])


# In[101]:


t = np.arange(256)
sp = np.fft.fft(vs[700])
efreq = np.fft.fftfreq(40)
plt.plot(efreq, sp.real, efreq, sp.imag)
plt.show()


# In[106]:


efreq


# In[82]:





# In[143]:


omegas[6000]/1e9

