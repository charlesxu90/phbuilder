#!/bin/python3

input  = open('MD_occ.pdb', 'r')
output = open('MD_occ_chain.pdb', 'w+')

chain  = ['A', 'B', 'C', 'D', 'E']
chainidx = 0
prev = 0

for line in input.readlines():

    type = line[0:6].strip()

    if type == 'MODEL':
        output.write(line)
        chainidx = 0
        prev = 0

    elif type in ['TITLE', 'REMARK', 'CRYST1', 'TER', 'ENDMDL']:
        output.write(line)

    elif type == 'ATOM':
        resname = line[17:21].strip()
        if resname in ['NA', 'CL', 'BUF', 'SOL', 'POPC']:
            output.write(line)
            continue

        output.write(line[0:21] + chain[chainidx] + line[22:])

        resid = int(line[22:27])
        if prev > resid:
            chainidx += 1
        prev = resid

input.close()
output.close()
