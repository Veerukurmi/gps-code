#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
from qutip import*
import openpyxl
from openpyxl import Workbook, load_workbook
import math as m
from qutip.piqs import *


# In[2]:


def Energy_Coherence(m,JK,gamma,delta,alpha,w0,w1):
    id2 = qeye(2)
    operators = [id2] * (m-1) 
    z=[tensor(*operators[:i], sigmaz(), *operators[i:]) for i in range(m)] 
    x=[tensor(*operators[:i], sigmax(), *operators[i:]) for i in range(m)] 
    y=[tensor(*operators[:i], sigmay(), *operators[i:]) for i in range(m)]  
    xx = []
    for i in range(len(x)-1):
        product =(1+gamma)* x[i]*x[i+1]+(1-gamma)*y[i]*y[i+1]+delta*z[i]*z[i+1]  #Nearest neighbour interaction
        xx.append(product)
    for i in range(len(x)-2):
        product =alpha*((1+gamma)*x[i]*x[i+2]+(1-gamma)*y[i]*y[i+2]+delta*z[i]*z[i+2]) #Next nearest neighbour interaction
        xx.append(product)
    H1=JK*(sum(xx)) #Interaction hamiltonian
    Hc=w1*(sum(x)) #Charging Hamiltonian
    HQ=w0*sum(z)  #Hamiltonian of the spins
    times = np.linspace(0,1.57,300) #times
    fS1 = [basis(2,1)]*m  
    psi1 = tensor(basis(2,1), *fS1[1:]) #state
    result = mesolve(Hc+H1, psi1, times,[],[HQ])
    E0=expect(HQ,psi1)
    Energy=result.expect[0]-E0
    result1 = mesolve(Hc+H1, psi1, times,[],[])
    Coherence=[]
    for i in range(len(result1.states)):
        rho1=result1.states[i]*result1.states[i].dag()
        diagonal_elements = np.diagonal(rho1)
        off_diagonal_elements = np.abs(np.array(rho1) - np.diag(diagonal_elements))
        coherence = np.sum(off_diagonal_elements)
        Coherence.append(coherence)
    return (Energy,Coherence)


# In[3]:


t=np.linspace(0,1.57,300)
Energy,Coherence= Energy_Coherence(10,0,0,0,0,1,1)


# In[20]:


Energy1,Coherence1= Energy_Coherence(10,1,0.2,1,1,1,1)
Energy2,Coherence2= Energy_Coherence(10,1,0.2,1,0,1,1)


# In[5]:


def Ins_pow(Energy,time):
    energy = np.array(Energy)
    time = np.array(t)
    delta_energy = np.diff(energy)
    delta_time = np.diff(time)
    IP= delta_energy / delta_time
    return IP


# In[7]:


P=Ins_pow(Energy,t)
np.max(P)


# In[13]:


plt.plot(t[1:],(1/19.99968504771061)*np.array(P),'g',linewidth=2)
# plt.plot(t[1:],(1/14.492105236990996)*np.array(P1),'k',linewidth=2)
# plt.plot(t[1:],(1/14.492105236990996)*np.array(P2),'r',linewidth=2)
plt.xlabel(r"$t/t_{min}$",fontsize=15)
plt.ylabel(r"$\mathcal{P}(t)/\mathcal{P}_{max}$",fontsize=15)
plt.xticks(np.linspace(0,1.57,6),[r'$0$',r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$',r'$1$'])


# In[35]:


plt.plot(t,Energy,'g',linewidth=2)
plt.xlabel(r"$t/t_{min}$",fontsize=15)
plt.ylabel(r"$\mathcal{E(t)}/\mathcal{E}_{max}$",fontsize=15)
plt.xticks(np.linspace(0,1.57,6),[r'$0$',r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$',r'$1$'])


# In[9]:


P=Energy/t
np.max(P[1:])


# In[16]:


plt.plot(t,(1/14.492105236990996)*np.array(P),'g',linewidth=2)
plt.xlabel(r"$t/t_{min}$",fontsize=15)
plt.ylabel(r"$\mathcal{P(t)}/\mathcal{P}_{max}$",fontsize=15)
plt.xticks(np.linspace(0,1.57,6),[r'$0$',r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$',r'$1$'])


# In[9]:


plt.plot(t,(1/1023)*np.array(Coherence),'g',linewidth=2)
# plt.plot(t,(1/1023)*np.array(Coherence1),'k',linewidth=2)
# plt.plot(t,(1/1023)*np.array(Coherence2),'r',linewidth=2)
plt.xlabel(r"$t/t_{min}$",fontsize=15)
plt.ylabel(r"$\mathcal{C(t)}$",fontsize=15)
plt.xticks(np.linspace(0,1.57,6),[r'$0$',r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$',r'$1$'])


# ## Varius results were comapered with the santos paper result. I found them compatible with my code.

# In[153]:


Pavg=[0.2333258685822128,0.23494675808524654,0.2687471118411693,0.3700904344371993,0.5141296743226476,0.6214810496866001,0.6369422594685555]
alpha=[-1,-0.8,-0.4,0,0.4,0.8,1]
plt.plot(alpha,Pavg,'b--')


# In[155]:


# Create a list starting from 0 to 1
my_list = np.linspace(0, 1, num=20)  # Adjust 'num' for the desired number of elements

# Print the list
print(my_list)


# ## Here I am writing codes for XYZ ordered interaction.

# ### Power for various anistropic parameters

# In[16]:


Energy1= Energy_Coherence(10,-1,0,1,0,1,1)
Energy2= Energy_Coherence(10,-1,0,1,0.5,1,1)
Energy3= Energy_Coherence(10,-1,0,1,0.8,1,1)


# In[12]:


t=np.linspace(0,1.57,300)


# In[17]:


Power=Energy1/t
np.max(Power[1:])


# In[14]:


# alpha=0
P0=[14.492095124830362,14.488865513140977,14.700711,15.335109489440036,16.365478585126798,17.74943467,19.496182781837284,21.615560580862784,24.05157075,26.720270629977964,29.55246548]
# alpha=0.5
P05=[14.49209858274007,14.669754944816576,15.38779259,16.61711700598563,18.309784520226323,20.4707972,23.020414006757733,25.854792501264306,28.88637329,32.054043984229644,35.31557854]
#alpha=0.8
P08=[14.4921027196962,14.800604281582139,15.86835431,17.46523726268363,19.582420246045682,22.15877434,25.061369756685636,28.21608376269163,31.57968393,35.10932923753174,38.76933267]
# Gamma
Gamma=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]


# In[15]:


plt.plot(Gamma,(1/14.492105236990996)*np.array(P0),'r-o',label=r"$\alpha=0$")
plt.plot(Gamma,(1/14.492105236990996)*np.array(P05),'k-o',label=r"$\alpha=0.5$")
plt.plot(Gamma,(1/14.492105236990996)*np.array(P08),'b-o',label=r"$\alpha=0.8$")
plt.xlabel(r"$\gamma$",fontsize=15)
plt.ylabel(r"$\overline{\mathcal{P}}_{max}/\mathcal{P}^{||}_{max}$", fontsize=15)
plt.legend()


# ### Varying N

# In[85]:


Energy1= Energy_Coherence(4,1,0.5,1,0.5,1,1)
Energy2= Energy_Coherence(5,1,0.5,1,0.5,1,1)
Energy3= Energy_Coherence(6,1,0.5,1,0.5,1,1)
Energy4= Energy_Coherence(7,1,0.5,1,0.5,1,1)
Energy5= Energy_Coherence(8,1,0.5,1,0.5,1,1)
Energy6= Energy_Coherence(9,1,0.5,1,0.5,1,1)
Energy7= Energy_Coherence(10,1,0.5,1,0.5,1,1)


# In[86]:


Energy8= Energy_Coherence(11,1,0.5,1,0.5,1,1)
Energy9= Energy_Coherence(12,1,0.5,1,0.5,1,1)
Energy10= Energy_Coherence(13,1,0.5,1,0.5,1,1)
Energy11= Energy_Coherence(14,1,0.5,1,0.5,1,1)
Energy12= Energy_Coherence(15,1,0.5,1,0.5,1,1)
Energy13= Energy_Coherence(16,1,0.5,1,0.5,1,1)
Energy14= Energy_Coherence(17,1,0.5,1,0.5,1,1)


# In[93]:


Power=Energy7/t
np.max(Power[1:])


# In[94]:


# When gamma =0
P00=[5.7968470042374145,7.246062248715683,8.695276296548679,10.144488199111516,11.593696860764364,13.042900931465924,14.49209858274007,15.941287241015763,17.39046343942144,18.83962235772177,20.2887572794221,21.737858797481014,23.18691368455708]
# When gamma=0.5
P005=[7.911113450864654,10.015496612374925,12.133618557589825,14.207634940194492,16.285928664778623,18.375969380637954,20.47079719646396]
N=[4,5,6,7,8,9,10]
plt.plot(N,P005)


# ### Varying J/h

# In[174]:


Energy11= Energy_Coherence(10,0.1,0.5,1,0,1,1)
Energy21= Energy_Coherence(10,0.1,0.5,1,0.5,1,1)
Energy31= Energy_Coherence(10,0.1,0.5,1,0.8,1,1)


# In[175]:


Power=Energy11/t
np.max(Power[1:])


# In[3]:


# alpha=0
P00=[14.567547712453607,14.91047324637575,15.360078871613863,15.88531417988523,16.465627232375027,17.089265858490787,17.74943466935851]
NP00=sorted(P00, reverse=True)+[14.369799743399836,14.339336908686018,14.428801765146476,14.492105236990996]
# alpha=0.5
P005=[14.47250900777186,14.537899266760357,14.831245276937388,15.343988232084333,16.00979559092507,16.776988884602474,17.618224077867435,18.51910784106935,19.471782434460774,20.47079719646396,14.492105236990996]
NP005=sorted(P005, reverse=True)
# alpha=0.8
P008=[14.490556521168038,14.665001417692263,15.140655515648103,15.842935715481648,16.689323328220954,17.637920311443388,18.66905866108984,19.772183201749264,20.938161150132952,22.158774344669077,14.492105236990996]
NP008=sorted(P008, reverse=True)
NN=[-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0]


# In[4]:


plt.plot(NN,(1/14.492105236990996)*np.array(NP00),'r-o',label=r"$\alpha=0$")
plt.plot(NN,(1/14.492105236990996)*np.array(NP005),'k-o',label=r"$\alpha=0.5$")
plt.plot(NN,(1/14.492105236990996)*np.array(NP008),'b-o',label=r"$\alpha=0.8$")
plt.xlabel(r"$J/\Omega$",fontsize=15)
plt.ylabel(r"$\overline{\mathcal{P}}_{max}/\mathcal{P}^{||}_{max}$", fontsize=15)
plt.legend()


# In[ ]:




