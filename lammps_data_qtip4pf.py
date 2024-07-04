import sys
import numpy as np
from ase.io import read

def get_lammps_data_qtip4pf(initfile):
    """
    Reads an input file (assumes OHH OHH ...) and
    returns a LAMMPS data file.
    """

    atoms = read(initfile)

    natoms = atoms.get_global_number_of_atoms()
    nbonds =  int(natoms / 3 * 2)
    nangles = int(natoms / 3)

    h = np.asarray(atoms.get_cell()).T / 0.529177
    q = np.asarray(atoms.get_positions()) / 0.529177
    s = atoms.get_chemical_symbols()
    s_i = [1 if tmp == "O" else 2 for tmp in s]

    qO = -1.1128
    qH = -qO / 2
    c_i = [qO if tmp == "O" else qH for tmp in s]

    print ("LAMMPS Description")
    print ("")
    print ("%10d atoms" % (natoms))
    print ("%10d bonds" % (nbonds))
    print ("%10d angles" % (nangles))
    print ("")
    print ("%10d atom types" % (2))
    print ("%10d bond types" % (1))
    print ("%10d angle types" % (1))
    print ("")
    print ("%10.5e %10.5e xlo xhi" % (0, h[0,0]))
    print ("%10.5e %10.5e ylo yhi" % (0, h[1,1]))
    print ("%10.5e %10.5e zlo zhi" % (0, h[2,2]))
    print ("%10.5e %10.5e %10.5e xy xz yz" % (h[1,0], h[2,0], h[2,1]))
    print ("")
    print ("Masses")
    print ("")
    print ("%10d %10.5e" %(1, 15.9994))
    print ("%10d %10.5e" %(2, 1.0080))
    print ("")
    print ("Bond Coeffs")
    print ("")
    print ("1    1.78    0.2708585 -0.327738785 0.231328959")
    print ("")
    print ("Angle Coeffs")
    print ("")
    print ("1    0.0700  107.400000")
    print ("")
    print ("Atoms")
    print ("")
    for i in range(natoms):
        # index_atom index_molecule index_atom_type charge x y z
        print ("%10d %10d %10d %10.5e %10.5e %10.5e %10.5e" % (i + 1 , i // 3 + 1, s_i[i], c_i[i], q[i][0], q[i][1], q[i][2]) )
    print ("")
    print ("Bonds")
    print ("")
    j = 0
    for i in range(0, nbonds // 2, 1):
        j += 1
        print ("%10d %10d %10d %10d" % (j, 1, (i * 3) + 1, (i * 3) + 2))
        j += 1
        print ("%10d %10d %10d %10d" % (j, 1, (i * 3) + 1, (i * 3) + 3))
    print ("")
    print ("Angles")
    print ("")
    for i in range(nangles):
        print ("%10d %10d %10d %10d %10d" % (i + 1, 1, (i * 3) + 2, (i * 3) + 1, (i * 3) + 3))


args = sys.argv[1:]
get_lammps_data_qtip4pf(args[0])
