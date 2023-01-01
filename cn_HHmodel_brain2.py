from brian2 import *

start_scope()

area = 20000 * umetre ** 2
Cm = 1 * ufarad * cm ** -2 * area
g_l = 0.03 * msiemens * cm ** -2 * area
El = -54 * mV
EK = -77 * mV
ENa = 50 * mV
g_Na = 120 * msiemens * cm ** -2 * area
g_k = 36 * msiemens * cm ** -2 * area


eqs = Equations('''
dv/dt = (I - (g_Na * m ** 3 * h * (v - ENa)) - (g_k * n ** 4 * (v - EK)) - (g_l * (v - El))) / Cm : volt

dn/dt = alpha_n * (1 - n) - beta_n * n : 1
dm/dt = alpha_m * (1 - m) - beta_m * m : 1
dh/dt = alpha_h * (1 - h) - beta_h * h : 1

alpha_n = (0.01 * (v / mV + 50) / (1 - exp(-(v / mV + 50) / 10))) / ms : Hz
beta_n = (0.125 * exp(-(v / mV + 60) / 80)) / ms : Hz

alpha_m = (0.1 * (v / mV + 35) / (1 - exp(-(v / mV + 35) / 10))) / ms : Hz
beta_m = (4.0 * exp(-0.0556 * (v / mV + 60))) / ms : Hz

alpha_h = (0.07 * exp(-0.05 * (v / mV + 60))) / ms : Hz
beta_h = (1 / (1 + exp(-(0.1) * (v / mV + 30)))) / ms : Hz

I : amp
g_Na : siemens
g_k : siemens
g_l : siemens
''')

group = NeuronGroup(1, eqs, method='exponential_euler')

group.v = -65 * mV
group.g_Na = g_Na
group.g_k = g_k
group.g_l = g_l
group.n = 0.15
group.m = 0.05
group.h = 0.75

state = StateMonitor(group, True, record=True)

group.I = '0*nA'
run(20 * ms)
group.I = '1*nA'
run(30 * ms)
group.I = '0*nA'
run(20 * ms)

plt.rcParams['figure.figsize'] = [12, 12]

plt.subplot(3,1,1)
plt.plot(state.t/ms, state.v[0]/mV)
plt.xlabel('Time (ms)')
plt.ylabel('v (mV)')

plt.subplot(3,1,2)
plt.plot(state.t/ms, state.m[0], label='m', color='red')
plt.plot(state.t/ms, state.n[0], label='n', color='green')
plt.plot(state.t/ms, state.h[0], label='h', color='blue')
plt.xlabel('Time (ms)')
plt.ylabel('Gating variable(probability)')
plt.legend()

plt.subplot(3,1,3)
plt.plot(state.t/ms, state.I[0]/nA)
plt.xlabel('Time (ms)')
plt.ylabel('I (nA)')