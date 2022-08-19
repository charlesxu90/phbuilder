#!/bin/python3

import analysis
import os
import pickle
import time
import sys

# PRE STUFF

if sys.argv[1] == 'real':
    dt = 20         # only use 1/20 frames
    b  = 0          # Do not ignore any frames, load entire trajectory.
    # b  = 300000   # ignore first 300'000 steps = 300 ns
elif sys.argv[1] == 'test':
    os.chdir('test')
    dt = 1
    b  = 2000
else:
    print("arguments are \'real\' or \'test\'")
    exit()

class Timer:
    def __init__(self, description=''):
        self.description = description
        self.__time1 = time.time()

    def time(self):
        time2 = time.time() - self.__time1
        if self.description == '':
            print('{:.3f}'.format(time2))
        else:
            print('{} took {:.3f} s'.format(self.description, time2))

# PICKLE PART

if not os.path.isfile('data.pickle'):

    print('No data.pickle detected, will construct and dump GLIC opbject...')

    t = Timer('Constructing GLIC object')
    GLIC = analysis.GLICSims(['4HFI_4','4HFI_7','6ZGD_4','6ZGD_7'], ['01', '02', '03', '04'], dt=dt, b=b)
    t.time()

    t = Timer('Dumping GLIC object')
    with open('data.pickle', 'wb') as handle:
        pickle.dump(GLIC, handle, protocol=pickle.HIGHEST_PROTOCOL)
    t.time()

else:
    print('Loading GLIC object from pickle data to save time...')
    
    t = Timer('Loading GLIC object from data')
    with open('data.pickle', 'rb') as handle:
        GLIC = pickle.load(handle)
    t.time()

# ACTUAL ANALYSIS CODE

t = Timer('Making plots')

GLIC.histograms(b=15000) # good
GLIC.convergence()
# GLIC.convergence_old(window=5000)
GLIC.histidine(b=15000)
GLIC.hisheatmap(b=15000)
GLIC.doFinalPlots()

t.time()
