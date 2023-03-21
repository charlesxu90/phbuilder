#!/usr/bin/env python3

import numpy as np

import MDAnalysis as mda
from MDAnalysis.analysis import distances

import argparse
import sys

def find_waters_to_remove(
    universe,
    cutoff,
    method,
    nNeighbours,
    minNWater,
    globalRadius,
    minNWaterAtoms
):
    if method == 'global':
        return global_waters_around_waters(
            universe,
            cutoff = cutoff,
            globalRadius = globalRadius,
            nmin = minNWaterAtoms
        )
    elif method == 'neighb':
        return neighbour_water_around_waters(
            universe,
            cutoff = cutoff,
            nn = nNeighbours,
            nmin = minNWater
        )
    else:
        raise ValueError("Wrong method provided")

def global_waters_around_waters(u, cutoff = 5, globalRadius = 10, nmin = 150):
    burried_resids = []
    resids = u.select_atoms(
        f"(resname SOL) and around {cutoff} (not resname SOL)"
    ).residues
    for r in resids:
        ids = r.atoms.ids
        nwateratoms = u.select_atoms(
            f"(resname SOL) and around {globalRadius} "
            + f"(id {ids[0]} {ids[1]} {ids[2]} and resid {r.resid})"
        ).n_atoms
        if nwateratoms <= nmin:
            burried_resids.append(r)
    return burried_resids

def neighbour_water_around_waters(u, cutoff = 5, nn = 20, nmin = 3):
    burried_resids = []
    residues = u.select_atoms(
        f"resname SOL and around {cutoff} (not resname SOL)"
    ).residues
    for r in residues:
        ids = r.atoms.ids
        resids = u.select_atoms(
            f"(around 5 (id {ids[0]} {ids[1]} {ids[2]} "
            + f"and resid {r.resid})) "
            + f"and (not (id {ids[0]} {ids[1]} {ids[2]} "
            + f"and resid {r.resid}))"
        ).residues
        dists = []
        for r1 in resids:
            dists.append(
                np.min(
                    distances.distance_array(
                        r.atoms.positions,
                        r1.atoms.positions
                    )
                )
            )
        residssorted = [x for _, x in sorted(zip(dists, resids),
            key=lambda pair: pair[0])]
        nwater = sum(1 for i in residssorted[:nn] if i.resname == "SOL")
        if nwater <= nmin:
            burried_resids.append(r)
    return burried_resids

def remove_from_structure_file(universe, waters_to_remove):
    if len(waters_to_remove) == 0:
        print("No waters have to be removed.")
        u.select_atoms("all").write("solvated_updated.gro")
        return 0
    ids = waters_to_remove[0].atoms.ids
    selection = (
        "not (resname SOL and (resid "
        + f"{waters_to_remove[0].resid} and id {ids[0]} {ids[1]} {ids[2]})"
    )
    for r in waters_to_remove[1:]:
        ids = r.atoms.ids
        selection += f" or (resid {r.resid} and id {ids[0]} {ids[1]} {ids[2]})"
    selection += ")"
    print(selection)
    u.select_atoms(selection).write("solvated_updated.gro")

def modify_topology(waters_to_remove):
    with open("topol.top", "r") as f, open("topol_updated.top", "w") as f1:
        for line in f:
            if line.strip().startswith("SOL"):
                nsol = int(line.split()[1])
                f1.write(f"SOL      {nsol-len(waters_to_remove)}\n")
            else:
                f1.write(line)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='This script will search water molecules that '
        + 'are accidentally placed inside the protein by gmx solvate.'
        + ' There are two algorithms available: (i) for each water '
        + 'molecule in the proximity of protein the number of waters'
        + ' within N nearest neighbours is computed (neighb method)'
        + ', if the number of water neighbours is smaller than fixed'
        + ' value M, the water molecule is considered misplaced and is'
        + ' deleted; (ii) for each water molecule in the proximity of '
        + 'protein the number of water atoms in the sphere of radius R '
        + 'is computed (global method), if number of water atoms is '
        + 'smaller than fixed value K, the water molecule is deleted.'
        + 'We provide the default values for N, M, R, and K, but the '
        + 'user can redefine it using the corresponding flags. Also, '
        + 'user can define the cutoff for the initial search of waters'
        + ' in the proximity of protein.'
    )

    parser.add_argument('-f', '--file', required=True, type=str,
            help='Input solvated structure (.gro/.pdb)')
    parser.add_argument('-p', '--top', required=True, type = str,
            help = 'Input topology file')
    parser.add_argument('-c', '--cutoff', required=False, type = float,
            default = 5,
            help = 'Cutoff (in angstroms) for the initial water search'
            + ' in the proximity of protein. Default is 5')
    parser.add_argument('-m', '--method', required=False, type = str,
            help = 'Method to use for searching misplaced waters. '
            + 'Default is global.',
            choices = ['global', 'neighb'], default = 'global')
    parser.add_argument('-N', '--nn', required=False, type = int,
            default = 20, help = 'Number of nearest neighbours.'
            + ' Default is 20.')
    parser.add_argument('-M', '--minn', required=False, type = int,
            default = 3, help = 'Maximum number of waters within '
            + 'nearest neighbours to remove water. Default is 3.')
    parser.add_argument('-r', '--gradius', required=False, type = float,
            default = 10, help = 'Sphere radius around water molecule '
            + 'within which the number of waters is calculated. Default'
            + ' is 10 angstroms.')
    parser.add_argument('-k', '--minatoms', required=False, type = int,
            default = 120, help ='Maximum number of water atoms within '
            + 'radius R of water molecules to remove water. Default is 120.')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    solvated = args.file
    u = mda.Universe(solvated)
    topol = args.top
    cutoff = args.cutoff
    method = args.method
    neigh = args.nn
    minneigh = args.minn
    globalr = args.gradius
    minatoms = args.minatoms

    waters_to_remove = find_waters_to_remove(
        u,
        cutoff,
        method,
        neigh,
        minneigh,
        globalr,
        minatoms
    )

    remove_from_structure_file(u, waters_to_remove)
