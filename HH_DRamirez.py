#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 14:48:43 2017

@author: davidramirez
"""

import scipy as sp
import pylab as plt
from scipy.integrate import odeint

# HH

C_m  =   1.0 
g_Na = 120.0
g_K  =  36.0
g_L  =   0.3
E_Na =  78.5
E_K  = -111.1
E_L  = -49.0

def alpha_m(V): 
    return -0.1*(V+35.0)/(sp.exp(-(V+35.0) / 10.0) - 1)

def beta_m(V):  
    return 4.0*sp.exp(-(V+60.0) / 18.0)

def alpha_h(V): 
    return 0.07*sp.exp(-(V+60.0) / 20.0)

def beta_h(V):  
    return 1.0/(1.0 + sp.exp(-(V+30.0) / 10.0))

def alpha_n(V): 
    return 0.01*(V+50.0)/(1.0 - sp.exp(-(V+50.0) / 10.0))

def beta_n(V):  
    return 0.125*sp.exp(-(V+60) / 80.0)


def I_Na(V,m,h):return g_Na * m**3 * h * (V - E_Na)

def I_K(V, n):  return g_K  * n**4     * (V - E_K)

def I_L(V):     return g_L             * (V - E_L)

# electrode current
def I_in(t):
    I_in = 10-10*(t>300)+40*(t>500)-40*(t>800)
    return I_in

# time
t = sp.arange(0.0, 1000.0, 0.1)

def ddt(dudes, t):
    V, m, h, n = dudes
    
    dVdt = (-1/C_m)*( I_Na(V, m, h) + I_K(V, n) + I_L(V) - I_in(t))
    dmdt = alpha_m(V)*(1.0-m) - beta_m(V)*m
    dhdt = alpha_h(V)*(1.0-h) - beta_h(V)*h
    dndt = alpha_n(V)*(1.0-n) - beta_n(V)*n
    return dVdt, dmdt, dhdt, dndt
    
X = odeint(ddt, [-83, 0.053, 0.6, 0.32], t)
V = X[:,0]
m = X[:,1]
h = X[:,2]
n = X[:,3]
ina = I_Na(V,m,h)
ik = I_K(V, n)
il = I_L(V)

plt.figure()

plt.title('Action potential due to different input current (Part A)')
plt.plot(t, V,'k')
plt.xlabel('time(ms)')
plt.ylabel('V (mV)')

#for a single action potential

I_in2 = 10
t2 = sp.arange(0.0, 20.0, 0.1)


def d2dt(dudes, t2):
    V, m, h, n = dudes
    
    dVdt = (-1/C_m)*( I_Na(V, m, h) + I_K(V, 5*n) + I_L(V) - I_in2)
    dmdt = alpha_m(V)*(1.0-m) - beta_m(V)*m
    dhdt = alpha_h(V)*(1.0-h) - beta_h(V)*h
    dndt = alpha_n(V)*(1.0-n) - beta_n(V)*n
    return dVdt, dmdt, dhdt, dndt

Y = odeint(ddt, [-83, 0,1,0],t2)
V2 = Y[:,0]
m2 = Y[:,1]
h2 = Y[:,2]
n2 = Y[:,3]
ina2 = I_Na(V,m,h)
ik2 = I_K(V, n)
il2 = I_L(V)

Z = odeint(d2dt, [-83, 0.053,0.6,0.32],t2)
V3 = Z[:,0]
m3 = Z[:,1]
h3 = Z[:,2]
n3 = Z[:,3]
ina3 = I_Na(V,m,h)
ik3 = I_K(V, n)
il3 = I_L(V)

plt.figure()
plt.title('Action potential (part B)')
plt.plot(t2,V2)
plt.ylabel('V (mV)')
plt.xlabel('time(ms)')

plt.figure()
plt.title('Gating probability(Part B)')
plt.plot(t2,m2,'r',label='$m$')
plt.plot(t2,n2,'g',label='$n$')
plt.plot(t2,h2,'b',label='$h$')
plt.ylabel('Prob')
plt.xlabel('time(ms)')
plt.legend()

plt.figure()
plt.title('Action potential with 5x n (part C)')
plt.plot(t2,V3)
plt.ylabel('V (mV)')
plt.xlabel('time(ms)')

plt.figure()
plt.title('Gating probability with 5x n (Part C)')
plt.plot(t2,m3,'r',label='$m$')
plt.plot(t2,n3,'g',label='$n$')
plt.plot(t2,h3,'b',label='$h$')
plt.ylabel('Prob')
plt.xlabel('time(ms)')
plt.legend()

#Question 1: Using Nernst potential equation:
    #E_K = -111.1 mV
    #E_Na = 78.5 mV
    #E_Cl = 45.6 mV
    
#Question 2: Using GHK equation:
    # Vm = -82.8 mV

#Question 3: [K]out increased 5 times:
    # Vm = -59.6 mV
    
#Question 4: Run code:
    #d: Increasing the gating variable n for potassium allows the channels
    # to remain open and then the membrane voltage will depend on mostly on K
    #since more potassium current will flow and overdrive the sodium potential
