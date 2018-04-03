# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 23:35:29 2018

@author: Chris
"""
from brian2 import *
tau_pre = tau_post = 20*ms
A_pre = 0.1
A_post = -A_pre*1.05
delta_t = linspace(-50, 50, 100)*ms
W = where(delta_t<0, A_pre*exp(delta_t/tau_pre), A_post*exp(-delta_t/tau_post))
plot(delta_t/ms, W)
xlabel(r'$\Delta t$ (ms)')
ylabel('$\Delta W$')
title('Connection Weighting Increment vs Time Difference')
ylim(-A_post, A_post)
axhline(0, ls='-', c='k')