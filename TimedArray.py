# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 21:13:23 2018

@author: Chris
"""

from brian2 import *
import matplotlib.pyplot as plt
from mido import Message, MidiFile, MidiTrack
import numpy as np
import mido

plt.close(0)
plt.close(1)
plt.close(2)
plt.close(3)

start_scope()


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

#visualise_connectivity(S)


# Initialise arrays for data storage
indexes = []
input_times = []
j = 0
test = 5

# Read in MIDI file times of notes on each track
mid = mido.MidiFile('John_Denver_-_Take_Me_Home_Country_Roads.mid')

MIDI_data = open('MIDI_data.txt', 'w')
for track in mid.tracks:   
    for msg in track:       
        MIDI_data.write("{0} \n".format(msg))   
    
for track in mid.tracks:
    for msg in track:
        if msg.type =="set_tempo":
            tempo = msg.tempo

slope_times = []

#tracks = [mid.tracks[4], mid.tracks[8]]
for track in mid.tracks:
    time = 0
    times = []
    #track = mid.tracks[4]
    for msg in track:
        if msg.type == "set_tempo":
            tempo = msg.tempo
        time+=(mido.tick2second(msg.time, mid.ticks_per_beat, tempo))
    
        if msg.type == "note_on":
            if msg.time == 0:
                time+=(mido.tick2second(1, mid.ticks_per_beat, tempo))
            input_times.append(time)
            times.append(time)
            indexes.append(j)       # Neuron indices
    #prev_time = msg.time
    slope_times.append(times)
    j+=1

# Initialise variables for use in the simulation
N = len(mid.tracks)
indices = array(indexes)
times = array(input_times)*second
taupre = 200*ms
taupost = taupre
dApre = 1
dApost = -dApre * taupre / taupost * 1.05
gmax = pi
#
## Represent the midi notes as spike times 
#inp = SpikeGeneratorGroup(N, indices, times)

runtime = 1*second
slopeArray = np.zeros(( N, int(runtime/100e-6)))

for neuron_index in range(N):
    timeprev = 0
    timedcounter = 0
    count = 0
    for j in range(len(slopeArray[neuron_index])):
        try:
            if timedcounter>slope_times[neuron_index][count]:
                timeprev = slope_times[neuron_index][count]
                count+=1
        except:
            continue
        slopeArray[neuron_index][j] = (((2*pi)/(slope_times[neuron_index][count]-timeprev)))
        timedcounter+=100e-6

#print(slopeArray)

slopeArrayTrans = np.transpose(slopeArray)
ta = TimedArray(slopeArrayTrans * second, dt=100*usecond)
#print(ta(100*ms,0))

# Represent the neurons as a set of equations
eqs_neurons = '''
dv/dt = ta(t,i) * (1/second)*(1/second) +  (c / N) * (1/second)  : 1
c : 1
'''

## Create a neuron group
G = NeuronGroup(N, eqs_neurons, threshold='v>2*pi', reset='v = 0', method='euler')
G.c = 0

# v_post - v_pre ???
recurrent = Synapses(G, G, '''w : 1
                dApre/dt = -Apre / taupre : 1 (clock-driven)
                dApost/dt = -Apost / taupost : 1 (clock-driven)''', 
                on_pre='''c_post = w*sin((v_pre-v_post))
                    Apre += dApre
                    w = clip(w + Apost, 0, gmax)''',
                on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''')


recurrent.connect('i!=j')
recurrent.w = 0
spikemon = SpikeMonitor(G)

trace = StateMonitor(G, ['v', 'c'], record=True)
trace0 = StateMonitor(recurrent, 'w', record=True)



run(runtime)
figure(0)
plot(trace.t/ms, trace.v.T)
xlabel('Time (ms)')
ylabel('Phase (rads)')
title('Oscillating Neurons With Coupling')

#legend(["Neuron 0", "Neuron 1"])
#plot(trace1.t/ms, trace1.v.T)
#plot(trace2.t/ms, trace2.v.T)

#figure(1)
#plot(trace.t/ms, trace.c.T)
#legend(j for j in range(len(mid.tracks)))
##plot(trace1.t/ms, trace1.m.T)
##plot(trace2.t/ms, trace2.m.T)
#
#figure(2)
##plot(trace0.t/ms, trace0.t_before.T)
#plot(spikemon.t/ms, spikemon.i, '.k')
#xlabel('Time (ms)')
#ylabel('Neuron index')
#title('Neuron Spike Timings')
#
#figure(3)
#plot(trace0.t/ms, trace0.w.T)
#xlabel('Time (ms)')
#ylabel('W')
#title('Connection weighting varying over time')
#legend(['Connection 0->1','Connection 1->0'])

figure(4)
visualise_connectivity(recurrent)

# %% Convert to MIDI

track_times_dict = {}
j = 0

for i in spikemon.i:
    track_times_dict[i] = []
    
for i in spikemon.i:
    track_times_dict[i].append(spikemon.t[j]/second)
    j+=1

i = 0


OUTPUT_data = open('Output_data.txt', 'w')

for track in mid.tracks:
    try:
        x = track_times_dict[i]
    except:
        print(i)
        i+=1
        continue
    j = 0    
    prev = 0
    for msg in track:
        if msg.type == "note_on":
            try: 
                #OUTPUT_data.write("Before: {0} \n".format(msg))
                OUTPUT_data.write("Prev: {0}\n".format(prev))
                OUTPUT_data.write("Track: {0}\n".format(track_times_dict[i][j]))
                OUTPUT_data.write("{0} \n".format(int(mido.second2tick((track_times_dict[i][j]-prev), mid.ticks_per_beat, tempo))))
                msg.time = int(mido.second2tick((track_times_dict[i][j]-prev), mid.ticks_per_beat, tempo))
                if msg.time<0:
                    msg.time *= -1
                #OUTPUT_data.write("After: {0} \n".format(msg))
                #prev = mido.tick2second(msg.time, mid.ticks_per_beat, tempo)
                prev += mido.tick2second(msg.time, mid.ticks_per_beat, tempo)
                j+=1
            except:
                prev += mido.tick2second(msg.time, mid.ticks_per_beat, tempo)
                j+=1
                continue
    i+=1
            
mid.save('timed_array.mid')
