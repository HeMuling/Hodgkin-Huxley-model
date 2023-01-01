'''
Hodgkin-Huxley model built on numpy

parameter function adpated from https://github.com/swharden/pyHH

method to solve the differential equation: Euler method
'''

import numpy as np

def alpha_n(v):
    return 0.01*(v+50)/(1-np.exp(-(v+50)/10))

def beta_n(v):
    return 0.125*np.exp(-(v+60)/80)

def alpha_m(v):
    return 0.1*(v+35)/(1-np.exp(-(v+35)/10))

def beta_m(v):
    return 4.0*np.exp(-0.0556*(v+60))

def alpha_h(v):
    return 0.07*np.exp(-0.05*(v+60))

def beta_h(v):
    return 1/(1+np.exp(-(0.1)*(v+30)))

def dmdt(V_pre, m_pre):

    return alpha_m(V_pre)*(1-m_pre)-beta_m(V_pre)*m_pre

def dndt(V_pre, n_pre):

    return alpha_n(V_pre)*(1-n_pre)-beta_n(V_pre)*n_pre

def dhdt(V_pre, h_pre):

    return alpha_h(V_pre)*(1-h_pre)-beta_h(V_pre)*h_pre

'''
current:
g_: in mS
m/h/n: no unit, probability
V_pre: in mV

return current: in uA
'''
def I_Na(V_pre, m_pre, h_pre, E_Na, g_Na):

    return g_Na * m_pre ** 3 * h_pre * (V_pre - E_Na)

def I_K(V_pre, n_pre, E_K, g_K):

    return g_K * n_pre ** 4 * (V_pre - E_K)

def I_L(V_pre, E_L, g_L):

    return g_L * (V_pre - E_L)

'''
voltage:
V_pre: in mV
time: in mS
I_e: in uA
C: in uF

return voltage: in mV
'''
def V_mem(V_pre, INa_pre, IK_pre, IL_pre, time_step, I_e, C=1):

    return V_pre + time_step * (I_e - INa_pre - IK_pre - IL_pre) / C

def cn_HHmodel(time, time_step, current_time, I_e, V_init,
                     ENa=50, EK=-77, EL=-54, g_Na=120, g_K=36,
                     g_L=0.03, Cm=1, area = 0.0002):
    '''
    time: total time to stimulate in mS
    time_step: in mS
    current_time: [start_time, end_time], in mS

    I_e: current we inject, in nA

    V_init: the inital membrane potential, in mV
    ENa, EK, EL: reversal potential, in mV

    g_Na, g_K, g_L: in mS/cm^2
    Cm: in uF/cm^2

    area: in cm^2
    '''

    # convert mS/cm^2 to mS; uF/cm^2 to uF
    g_Na *= area
    g_K *= area
    g_L *= area
    Cm *= area

    # initialize time
    t = np.linspace(0, time, int(time / time_step))
    length = len(t)

    begin_time, end_time = current_time # get the start time and end time for current to have effect
    
    # initialize parameter array
    v = np.zeros(length)
    m = np.zeros(length)
    h = np.zeros(length)
    n = np.zeros(length)

    # give inital conditions
    v[0] = V_init
    m[0] = alpha_m(v[0]) / (alpha_m(v[0]) + beta_m(v[0]))
    n[0] = alpha_n(v[0]) / (alpha_n(v[0]) + beta_n(v[0]))
    h[0] = alpha_h(v[0]) / (alpha_h(v[0]) + beta_h(v[0]))

    # compute the membrane potential
    dt = time_step
    for i in range(1,len(t)):

        if begin_time / 10 * 10000 <= i <= end_time / 10 * 10000:
            I = I_e * 1e-03
        else:
            I = 0
        
        m[i] = m[i-1] + dt * dmdt(v[i-1], m[i-1])
        n[i] = n[i-1] + dt * dndt(v[i-1], n[i-1])
        h[i] = h[i-1] + dt * dhdt(v[i-1],h[i-1])

        INa = I_Na(v[i-1], m[i-1], h[i-1],ENa, g_Na)
        IK = I_K(v[i-1], n[i-1], EK, g_K)
        IL = I_L(v[i-1], EL, g_L)

        v[i] = V_mem(v[i-1], INa, IK, IL, dt, I, Cm)

    return v, m, h, n, t

'''
an example
'''


import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = [14, 14]

v, m, h, n, t = cn_HHmodel(40, 0.001, [10, 13], 10, -65)

plt.subplot(2,1,1)
plt.plot(t, v)
plt.title('voltage')
plt.xlabel('time (ms)')
plt.ylabel('voltage (mV)')

plt.subplot(2,1,2)
plt.plot(t, m, label='m')
plt.plot(t, h, label='h')
plt.plot(t, n, label='n')
plt.title('gating variables')
plt.xlabel('time (ms)')
plt.ylabel('gating variables (prabability)')
plt.legend()