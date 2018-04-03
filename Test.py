# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 17:50:52 2018

@author: Chris
"""

from brian2 import *
import pylab as plt
from mido import Message, MidiFile, MidiTrack
import mido
# IPython import get_ipython
#get_ipython().run_line_magic('matplotlib', 'inline')

plt.close()
plt.close()
plt.close()

start_scope()

#%% 

def visualise_connectivity(S):
    Ns = len(S.source)
    Nt = len(S.target)
    figure(figsize=(10, 4))
    subplot(121)
    plot(zeros(Ns), arange(Ns), 'ok', ms=10)
    plot(ones(Nt), arange(Nt), 'ok', ms=10)
    for i, j in zip(S.i, S.j):
        plot([0, 1], [i, j], '-k')
    xticks([0, 1], ['Source', 'Target'])
    ylabel('Neuron index')
    xlim(-0.1, 1.1)
    ylim(-1, max(Ns, Nt))
    subplot(122)
    plot(S.i, S.j, 'ok')
    xlim(-1, Ns)
    ylim(-1, Nt)
    xlabel('Source neuron index')
    ylabel('Target neuron index')


#%% 



n = 6

eqs = '''
dv/dt = (w)*(1/second) + (K / n * sin(vj - v))*(1/second)  : 1
w : 1
vj : 1
K : 1
'''
neurons = NeuronGroup(n, model=eqs, threshold='v > 1', reset='v = 0',
                      method='euler')
neurons.w = [1, 1.03, 2, 2.04, 3, 3.07]
neurons.K = 0.8
neurons.v = 0

S = Synapses(neurons, neurons, on_pre = 'vj = v + rand()*0.1')
S.connect(condition='i!=j')

spikemon = SpikeMonitor(neurons)

trace0 = StateMonitor(neurons, 'v', record=0)
trace1 = StateMonitor(neurons, 'v', record=1)
trace2 = StateMonitor(neurons, 'v', record=2)
trace3 = StateMonitor(neurons, 'v', record=3)
trace4 = StateMonitor(neurons, 'v', record=4)
trace5 = StateMonitor(neurons, 'v', record=5)


run(10000*ms)
subplot(211)
plot(spikemon.t/ms, spikemon.i, '.k')
xlabel('Time (ms)')
ylabel('Neuron index')
subplot(212)
plot(trace0.t/ms, trace0.v.T)
plot(trace1.t/ms, trace1.v.T)
plot(trace2.t/ms, trace2.v.T)
plot(trace3.t/ms, trace3.v.T)
plot(trace4.t/ms, trace4.v.T)
plot(trace5.t/ms, trace5.v.T)
xlabel('Time (ms)')
ylabel('v')
tight_layout()
show()

#%% Kuramoto Sonification

track_times_dict = {}
j = 0

for i in spikemon.i:
    track_times_dict[i] = []
    
for i in spikemon.i:
    track_times_dict[i].append(spikemon.t[j])
    j+=1

mid = MidiFile()

program_counter = 0
for i in track_times_dict:
    track = MidiTrack()
    track.append(Message('program_change', program=program_counter, time=0))
    temp = 0
    for time in track_times_dict[i]:
        track.append(Message('note_on', note=64+(program_counter*2), velocity=64, time=(int(mido.second2tick(time - temp, mid.ticks_per_beat, 500000)))))
        track.append(Message('note_off', note=64+(program_counter*2), velocity=127, time=64*1))
        temp = time 
    mid.tracks.append(track)
    program_counter+=1

mid.save('new_song.mid')