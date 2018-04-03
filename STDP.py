# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 19:21:26 2018

@author: Chris
"""

plt.close()
plt.close()

from brian2 import *
from mido import Message, MidiFile, MidiTrack
import mido

start_scope()

# Determines how quickly neurons change value (do not have to be used)
taupre = 20*ms
taupost = taupre

# Number of oscillators
n = 3

# Determines maximum increase in K value
gmax = 2

# Initial conditions for weight change???
dApre = 0.7
dApost = -dApre * taupre / taupost * 1.05
#dApost *= gmax
#dApre *= gmax


# Equations to represent each neuron
# Each neuron has a phase (v), frequency (freq), coupling constant (K)
# vj is the phase variable of the target oscillator

# Should sine and coupling be implemented in the Synapse???
eqs_neurons = '''
dv/dt = freq * (1/second) + (c / n)*(1/second): 1
c : 1
freq : 1
'''


# Define a neuron group comprising of neurons represented by equations above
neurons = NeuronGroup(n, eqs_neurons, threshold='v>1', reset='v = 0',
                      method='euler')

neurons.freq = [2, 1.5, 2.3]


S = Synapses(neurons, neurons,
             '''w : 1
                dApre/dt = -Apre / taupre : 1 (clock-driven)
                dApost/dt = -Apost / taupost : 1 (clock-driven)''',
             on_pre='''c_post = w*sin((v_pre-v_post)*2*pi)
                    Apre += dApre
                    w = clip(w + Apost, 0, gmax)''',
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',
             )



S.connect(condition='i!=j')

neuron_mon = StateMonitor(neurons, ['v', 'c'], record=True)
synapse_mon = StateMonitor(S, ['w'], record = True)
spikemon = SpikeMonitor(neurons)

run(50*second, report='text')

subplot(311)
plot(S.v / gmax, '.k')
ylabel('Weight / gmax')
xlabel('Synapse index')
subplot(312)
hist(S.v / gmax, 20)
xlabel('Weight / gmax')
subplot(313)
plot(synapse_mon.t/second, synapse_mon.w.T/gmax)
xlabel('Time (s)')
ylabel('Weight / gmax')
tight_layout()
show()

figure()
#plot(neuron_mon.t/ms, neuron_mon.c.T)
plot(neuron_mon.t/ms, neuron_mon.v.T)


#%% Kuramoto Sonification

track_times_dict = {}
j = 0

for i in spikemon.i:
    track_times_dict[i] = []
    
for i in spikemon.i:
    track_times_dict[i].append(spikemon.t[j])
    j+=1

mid = MidiFile(type=1)

program_counter = 0

for i in track_times_dict:
    track = MidiTrack()
    track.append(Message('program_change', channel=9, program=program_counter, time=0))
    temp = 0
    for note_time in track_times_dict[i]:
        track.append(Message('note_on', channel=9, note=64+(program_counter*2), velocity=64, time=(int(mido.second2tick((note_time - temp - (mido.tick2second(64, mid.ticks_per_beat, 500000)*(second))), mid.ticks_per_beat, 500000)))))
        track.append(Message('note_off', channel=9, note=64+(program_counter*2), velocity=127, time=64))
        temp = note_time
    mid.tracks.append(track)
    program_counter+=1
    
mid.save('new_song.mid')