import os
import sys
import MDAnalysis

def lastLambdaVal(fname):
    """Retrieves the last lambda coordinate from a lambda file."""
    with open(fname) as file:
        for line in file:
            pass
        return float(line.split()[1])

base = sys.argv[1]
new  = sys.argv[2]

if not os.path.exists(base):
    os.mkdir(base)

# Extract and move the lambda-coordinate files.
os.system(f'gmx cphmd -s {base}.tpr -e {base}.edr -numplot 1')
os.system(f'mv cphmd-coord-*.xvg {base}')

# Compile a list of last lambda values for each lambda coordinate.
values = []
for num in range(1, len(os.listdir(base)) + 1):
    values.append(lastLambdaVal(f"{base}/cphmd-coord-{num}.xvg"))

# Compile a list of titratable residue names in the protein.
u = MDAnalysis.Universe('EM.gro')
acidics = list(u.select_atoms('resname ASPT GLUT HSPT').residues.resnames)
acidics.append('BUF')

# Create the strings to replace the current lines in upcoming .mdp with.
lidx = 0
for atomset in range(0, len(acidics)):

    if acidics[atomset] == 'HSPT':
        newCoord = str(values[lidx]) + ' ' + str(values[lidx + 1]) + ' ' + str(values[lidx + 2])
        lidx += 3

    if acidics[atomset] in ['ASPT', 'GLUT', 'BUF']:
        newCoord = str(values[lidx])
        lidx += 1

    print(atomset + 1, acidics[atomset], newCoord)

    # In the upcoming .mdp file, replace the line with the last coordinate.
    match   = f"lambda-dynamics-atom-set{atomset + 1}-initial-lambda"
    replace = match + f"               = {newCoord}"

    os.system(f'sed -i \'s/.*{match}.*/{replace}/\' {new}.mdp')
