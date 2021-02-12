#!/usr/bin/env python

import os
import collections
import argparse
import numpy
import spglib
import math

from . import crystal
from . import spg


def main():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
+-+-+-+-+-+-+-+-+-+-+-+-+
       C I F L I B
+-+-+-+-+-+-+-+-+-+-+-+-+
author : PearCandy
last update : 2021.2.12
contact y.nishi1980@gmail.com''')

    arg.add_argument("-c", "--cif", action="store", dest="cif")
    arg.add_argument("--poscar", action="store", dest="poscar")
    arg.add_argument("-p",
                     "--program",
                     action="store",
                     dest="program",
                     default="qe")
    arg.add_argument("-o", "--output", action="store", dest="fout")
    arg.add_argument("--ocif", action="store", dest="foutcif")
    arg.add_argument("--cif_calc", action="store", dest="cif_calc")
    arg.add_argument("-bz", "--bz_gnuplot", action="store", dest="bz_gnuplot")
    arg.add_argument("--bz_jmol", action="store", dest="bz_jmol")
    arg.add_argument("--nkpath", action="store", dest="nkpath", default=100)
    arg.add_argument("--calc", action="store", dest="calc")
    arg.add_argument("--tol", action="store", dest="tol", default="1.0e-5")
    args = arg.parse_args()
    cryst = crystal.Crystal()
    print("+-+-+-+-+-+-+-+-+-+-+-+-+")
    print("       C I F L I B       ")
    print("       ver.", cryst.version)
    print("    [usage] ciflib -h   ")
    print("")
    print("    author: PearCandy")
    print("  y.nishi1980@gmail.com")
    print("+-+-+-+-+-+-+-+-+-+-+-+-+")
    spgnum_cif = 0

    tol = float(args.tol)

    #cif exist or not?
    if args.cif != None:
        print("\n[load cif files]")
        print('cif:', args.cif)
        cryst.read_cif(args.cif)
        cryst.check_cif()

        cryst.construct_structure()

        print('\n[space group]')

        print(cryst.space_group_crystal_system)
        print(cryst.symmetry_Int_Tables_number)
        print(cryst.symmetry_space_group_name)

        print(cryst.cell.length_a, cryst.cell.length_b, cryst.cell.length_c)
        print(cryst.cell.angle_alpha, cryst.cell.angle_beta,
              cryst.cell.angle_gamma)

        spgnum_cif = cryst.symmetry_Int_Tables_number

        print('\n---- ciflib execution ----')
        b1, hallnum_ccell, hallnum_pcell = spg.check_by_ciflib(cryst, tol)
        print('ciflib>:', b1)
    else:
        print("\nPlease specify the CIF file\n")

    #poscar exist or not?
    if (args.poscar != None):
        print('poscar:', args.poscar)
        cryst.read_structure_vasp(args.poscar)
        cryst.reconstruct_structure()

        print('\n[calc]')
        b2 = spg.check_by_ciflib(cryst, tol)
        print('calc:', b2)

        cell, atoms = read_structure_vasp(args.poscar)

        print('\n[calc]')
        b3, spgnum_calc = calc_check(cell, spgnum_cif, tol)
        print('calc:', b3)

        if (b3 == True):
            print('structure is not changed.  orginal: %d  calc: %d' %
                  (spgnum_cif, spgnum_calc))
        else:
            print('warning: structure is changed.  orginal: %d  calc: %d' %
                  (spgnum_cif, spgnum_calc))

        if (args.cif_calc != None):
            cifdata = make_cifdata(cell, atoms, tol)
            print('output cif (calc):', args.cif_calc)
            write_cif_file(args.cif_calc, cifdata)

    #calc or not?
    if args.calc != None:
        print('calc:', args.calc)
        if args.program == "vasp":
            cell, atoms = read_structure_vasp(args.calc)
        elif args.program == "qe":
            cell, atoms = read_structure_qe(args.calc)
        elif args.program == "abinit":
            cell, atoms = read_structure_abinit(args.calc)

        if args.cif != None:
            print('[calc]')
            b3, spgnum_calc = calc_check(cell, spgnum_cif, tol)
            print('calc:', b3)

            if (b3 == True):
                print('structure is not changed.  orginal: %d  calc: %d' %
                      (spgnum_cif, spgnum_calc))
            else:
                print('warning: structure is changed.  orginal: %d  calc: %d' %
                      (spgnum_cif, spgnum_calc))

        if args.cif_calc != None:
            cifdata = make_cifdata(cell, atoms, tol)
            print('output cif (calc):', args.cif_calc)
            write_cif_file(args.cif_calc, cifdata)

    #nkpath exist or not?
    total_sampling = 100
    if args.nkpath != None:
        total_sampling = int(args.nkpath)

    #cif exist or not?
    if args.cif != None:
        if (args.program == 'vasp'):
            print('output: POSCAR, KPOINTS')
            cryst.write_structure_vasp('POSCAR')
            cryst.write_kpoints_vasp('KPOINTS', total_sampling)
            cryst.write_kpoints_vasp_line("KPOINTS.line", 10)
        elif (args.program == 'qe'):
            if (args.fout != None):
                print('output:', args.fout)
                cryst.write_structure_espresso(args.fout)
                cryst.write_kpoints_espresso(args.fout, total_sampling)
        elif (args.program == 'abinit'):
            if (args.fout != None):
                print('output:', args.fout)
                cryst.write_structure_abinit(args.fout)
                cryst.write_kpoints_abinit(args.fout, total_sampling)

        if (args.foutcif != None):
            print('output cif (update):', args.foutcif)
            cryst.write_cif(args.foutcif)

        if (args.bz_gnuplot != None):
            print('output bz (gnuplot):', args.bz_gnuplot)
            cryst.write_brillouinzone_gnuplot(args.bz_gnuplot)
        if (args.bz_jmol != None):
            print('output bz (jmol):', args.bz_jmol)
            cryst.write_brillouinzone_jmol(args.bz_jmol)


def calc_check(cell, spgnum, tol=1.0e-5):
    print('\n[primitive cell (calc)]')
    print_cell(cell)
    dataset = spglib.get_symmetry_dataset(cell, tol)
    spgnum_pcell = dataset['number']
    hallnum_pcell = dataset['hall_number']
    print_dataset(dataset)

    print('\n[conventional cell (calc pcell+standardize)]')
    scell = spglib.standardize_cell(cell)
    print_cell(scell)
    dataset = spglib.get_symmetry_dataset(scell, tol)
    spgnum_ccell = dataset['number']
    hallnum_ccell = dataset['hall_number']
    print_dataset(dataset)

    print('\n# check crystal structure...')
    print('%-25s %10s %10s' % ('', 'spg num', 'hall num'))
    print('%-25s %10d %10s' % ('cell orginal cif', spgnum, '-'))
    print('%-25s %10d %10d' % ('pcell calc', spgnum_pcell, hallnum_pcell))
    print('%-25s %10d %10d' %
          ('ccell calc pcell+std', spgnum_ccell, hallnum_ccell))
    b = True
    if (spgnum != spgnum_pcell):
        b = False
        print('warning: primitive cell... NG')
    if (spgnum != spgnum_ccell):
        b = False
        print('warning: conventional cell... NG')
    if (b == True):
        print('calc... OK')

    return b, spgnum_pcell


def read_structure_vasp(filename):
    input = open(filename, "r")
    input.readline()  # comment
    scale = float(input.readline())
    lattice = []
    for i in range(3):
        v0 = input.readline().strip().split()
        v = list(map(float, [v0[0], v0[1], v0[2]]))
        lattice.append([v[0] * scale, v[1] * scale, v[2] * scale])
    atoms = input.readline().strip().split()
    numatoms = list(map(int, input.readline().strip().split()))
    numtot = sum(numatoms)

    numbers = []
    for i in range(len(numatoms)):
        for j in range(numatoms[i]):
            numbers.append(i)
    input.readline()  # Direct assumed
    positions = []
    for i in range(numtot):
        v0 = input.readline().strip().split()
        v = list(map(float, [v0[0], v0[1], v0[2]]))
        positions.append(v)
    input.close()

    return (lattice, positions, numbers), atoms


def read_structure_qe(filename):
    f = open(filename, "r")
    lattice = []
    positions = []
    numbers = []
    atoms = []
    ele = {}
    for l in f:
        if (l.find('nat') >= 0):
            l0 = l.split('=')
            natom = int(l0[1])
        if (l.find('ntyp') >= 0):
            l0 = l.split('=')
            ntype = int(l0[1])
        if (l.find('CELL_PARAMETERS') >= 0):
            for i in range(3):
                c = f.next()
                l0 = c.split()
                v = list(map(float, l0))
                lattice.append([v[0], v[1], v[2]])
        if (l.find('ATOMIC_SPECIES') >= 0):
            for i in range(ntype):
                c = f.next()
                l0 = c.split()
                atoms.append(l0[0])
                ele[l0[0]] = i
        if (l.find('ATOMIC_POSITIONS') >= 0):
            for i in range(natom):
                c = f.next()
                l0 = c.split()
                v = list(map(float, [l0[1], l0[2], l0[3]]))
                positions.append(v)
                numbers.append(ele[l0[0]])

    return (lattice, positions, numbers), atoms


def read_structure_abinit(filename):
    f = open(filename, "r")
    lattice = []
    positions = []
    numbers = []
    atoms = []
    for l in f:
        if (l.find('acell') >= 0):
            l0 = l.split()
            acell = [float(l0[1]), float(l0[2]), float(l0[3])]
        if (l.find('rprim') >= 0):
            rprim = []
            for i in range(3):
                c = f.next()
                l0 = c.split()
                v = list(map(float, l0))
                rprim.append(v)
        if (l.find('ntypat') >= 0):
            l0 = l.split()
            ntype = int(l0[1])
        if (l.find('znucl') >= 0):
            l0 = l.split()
            znucl = []
            for i in range(ntype):
                znucl.append(l0[i + 1])
        if (l.find('natom') >= 0):
            l0 = l.split()
            natom = int(l0[1])
        if (l.find('typat') >= 0):
            l0 = l.split()
            for i in range(ntype):
                numbers.append(int(l0[i + 1]) - 1)
        if (l.find('xred') >= 0):
            for i in range(natom):
                c = f.next()
                l0 = c.split()
                v = list(map(float, l0))
                positions.append(v)
    f.close()

    lattice = numpy.multiply(rprim, acell)

    elements = {}
    cfgdir = str(os.path.dirname(__file__)) + "/config"
    f = open(cfgdir + "/elements.aw")
    for l in f:
        if "#" in l:
            continue
        l0 = l.split()
        atomno = l0[0]
        e = l0[1]
        elements[atomno] = e
    for a in znucl:
        atoms.append(elements[a])

    return (lattice, positions, numbers), atoms


def make_cifdata(cell, atoms, tol=1.0e-5):
    cifdata = []

    scell = spglib.standardize_cell(cell)
    print('*** scell ***')
    print_cell(scell)
    dataset = spglib.get_symmetry_dataset(scell, tol)

    wyckoffs = dataset['wyckoffs']
    equiv_atoms = dataset['equivalent_atoms']
    trans_matrix = dataset['transformation_matrix']
    origin_shift = dataset['origin_shift']
    hall_number = dataset['hall_number']
    hall_symbol = dataset['hall']
    sg_type = spglib.get_spacegroup_type(hall_number)
    print('trans_matrix')
    print(trans_matrix)
    print('origin_shift')
    print(origin_shift)

    lattice = scell[0]
    positions = scell[1]
    numbers = scell[2]
    cifdata.append(('data_', ''))
    cifdata.append(('', None))

    chemform = ''
    num = [0 for i in range(len(atoms))]
    for i, a in enumerate(atoms):
        for n in numbers:
            if (i == n):
                num[i] += 1
    if ('C' in atoms):
        chemform += 'C' + str(num[atoms.index('C')]) + ' '
    if ('H' in atoms):
        chemform += 'H' + str(num[atoms.index('H')]) + ' '
    for i, a in enumerate(atoms):
        if (a != 'C' and a != 'H'):
            chemform += a + str(num[i]) + ' '
    chemform = chemform.strip()
    print('chemical formula', chemform)

    cifdata.append(('_chemical_formula_sum', chemform))
    crystal_system = get_crystal_system(dataset['number'])
    cifdata.append(('_space_group_crystal_system', crystal_system))

    name_HM = sanitize_name_HM(sg_type['international_short'])
    cifdata.append(('_symmetry_space_group_name_H-M', name_HM))
    cifdata.append(('_symmetry_Int_Tables_number', dataset['number']))
    cifdata.append(('_symmetry_space_group_name_Hall', hall_symbol))
    cifdata.append(('', None))

    cifdata_symmetry_list = ('loop_symmetry', [])
    cifdata_symmetry_list[1].append(
        ['_symmetry_equiv_pos_site_id', '_symmetry_equiv_pos_as_xyz'])

    ind = 0
    for rot, trans in zip(dataset['rotations'], dataset['translations']):
        symbols = make_symop(rot, trans)
        ind += 1
        cifdata_symmetry_list[1].append([ind, symbols])

    cifdata.append(cifdata_symmetry_list)
    cifdata.append(("", None))

    a, b, c = celllength(lattice)
    alpha, beta, gamma = cellangle(lattice)
    cifdata.append(('_cell_length_a', a))
    cifdata.append(('_cell_length_b', b))
    cifdata.append(('_cell_length_c', c))
    cifdata.append(('_cell_angle_alpha', alpha))
    cifdata.append(('_cell_angle_beta', beta))
    cifdata.append(('_cell_angle_gamma', gamma))
    cifdata.append(('', None))

    # find inequivalent atoms and calculate their multiplicitis
    inequiv_atoms = collections.OrderedDict()
    for num in equiv_atoms:
        if num not in inequiv_atoms:
            inequiv_atoms[num] = 1
        else:
            inequiv_atoms[num] += 1

    print('space group number %s' % dataset['number'])
    print('space group %s' % dataset['international'])
    print('hall symbol %s' % dataset['hall'])

    print('name wyckoff  multiplicity  x y z')

    cifdata_atomlist = ('loop_atom', [])
    cifdata_atomlist[1].append([
        '_atom_site_label', '_atom_site_fract_x', '_atom_site_fract_y',
        '_atom_site_fract_z', '_atom_site_occupancy',
        '_atom_site_symmetry_multiplicity', '_atom_site_Wyckoff_symbol',
        '_atom_site_type_symbol'
    ])
    suffix_list = {}
    for atom in atoms:
        suffix_list[atom] = 1
    for num in inequiv_atoms:
        atom_type_number = numbers[num]
        atom_type = atoms[atom_type_number]
        atom_label = atom_type + str(suffix_list[atom_type])
        suffix_list[atom_type] += 1
        cifdata_atomlist[1].append([
            atom_label, positions[num][0], positions[num][1],
            positions[num][2], 1.0, inequiv_atoms[num], wyckoffs[num],
            atom_type
        ])
        print(atom_label, atom_type, wyckoffs[num], inequiv_atoms[num],
              positions[num][0], positions[num][1], positions[num][2])

    cifdata.append(cifdata_atomlist)
    cifdata.append(('', None))

    return cifdata


def write_cif_file(fname, cifdata):
    ''' write cif data to file
        fname: output filename
        cifdata:  cif data as a list of tuple (key, val)
                  key corresponds to a cif key.
                  if key is not loop_, print "key val"
                  if key is loop_, val[0] is a list of identifier in loop_ and val[1:] contains data.
    '''
    with open(fname, "w") as f:

        for key, val in cifdata:
            if key == 'data_':
                f.write('%s%s\n' % (key, val))
            elif key == "_chemical_formula_sum":
                f.write("%s '%s'\n" % (key, val))
            elif key == '_space_group_crystal_system':
                f.write('%s %s\n' % (key, val))
            elif key == '_symmetry_space_group_name_H-M':
                f.write("%s '%s'\n" % (key, val))
            elif key == '_symmetry_Int_Tables_number':
                f.write("%s %s\n" % (key, val))
            elif key == '_symmetry_space_group_name_Hall':
                f.write("%s '%s'\n" % (key, val))

            elif key == 'loop_symmetry':
                f.write('loop_\n')
                loop_keys = val[0]
                loop_vals = val[1:]
                for lk in loop_keys:
                    f.write('    %s\n' % lk)
                for lv in loop_vals:
                    f.write(" %d '%s'\n" % (lv[0], lv[1]))

            elif key.find('_cell_length') >= 0:
                f.write('%s %15.9f\n' % (key, val))
            elif key.find('_cell_angle') >= 0:
                f.write('%s %15.9f\n' % (key, val))

            elif key == 'loop_atom':
                f.write('loop_\n')
                loop_keys = val[0]
                loop_vals = val[1:]
                for lk in loop_keys:
                    f.write('    %s\n' % lk)
                for lv in loop_vals:
                    f.write(' %s %15.9f %15.9f %15.9f %f %s %s %s\n' % \
                                (lv[0],lv[1],lv[2],lv[3],lv[4],lv[5],lv[6],lv[7]))

            elif key == '':
                f.write('\n')

            else:
                f.write('%s %s' % (key, val))


def print_cell(cell):
    lattice, positions, numbers = cell
    print('lattice vectors')
    for vec in lattice:
        print('%15.9f %15.9f %15.9f' % (vec[0], vec[1], vec[2]))
    print('atomic postions and type')
    for coor, atype in zip(positions, numbers):
        print('%4d %15.9f %15.9f %15.9f' % (atype, coor[0], coor[1], coor[2]))


def print_dataset(dataset):
    print('number', dataset['number'])
    print('international', dataset['international'])
    print('hall_symbol', dataset['hall'])
    print('hall_number', dataset['hall_number'])
    print('main_axis_setting', dataset['choice'])
    print('pointgroup', dataset['pointgroup'])
    hall_number = dataset['hall_number']
    sg_type = spglib.get_spacegroup_type(hall_number)
    print_spgtype(sg_type)
    print('trans_matrix\n', dataset['transformation_matrix'])
    print('original_shift\n', dataset['origin_shift'])


def print_spgtype(sgtype):
    print('number', sgtype['number'])
    print('international_short', sgtype['international_short'])
    print('international_full', sgtype['international_full'])
    print('international', sgtype['international'])
    print('hall_sysmbol', sgtype['hall_symbol'])
    print('main_axis_setting', sgtype['choice'])
    print('pointgroup', sgtype['pointgroup_schoenflies'])


def sanitize_name_HM(name):
    work = []
    cur = ''
    length = len(name)
    i = 0
    while i < length:
        if name[i] not in '_/':
            if cur:
                work.append(cur)
            if name[i] == '-':
                cur = name[i:i + 2]
                i += 2
            else:
                cur = name[i]
                i += 1
        elif name[i] == '_':
            cur += name[i + 1]
            i += 2
        elif name[i] == '/':
            cur += name[i:i + 2]
            i += 2
    if cur:
        work.append(cur)
    return ' '.join(work)


def get_crystal_system(it_number):
    '''retrun crystal system from IT number'''
    if it_number <= 0:
        return 'error'
    elif it_number < 3:
        return 'triclinic'
    elif it_number < 16:
        return 'monoclinic'
    elif it_number < 75:
        return 'orthorhombic'
    elif it_number < 143:
        return 'tetragonal'
    elif it_number < 168:
        return 'trigonal'
    elif it_number < 195:
        return 'hexagonal'
    elif it_number < 231:
        return 'cubic'
    else:
        return 'error'


def make_symop(rot, trans):
    '''make symop symbols from rotation matrix and translation vectors
       rot: rotational matrix in numpy 3x3 array
       trans: translational vector
       return string 
    '''
    symbols = []
    for i in range(3):
        symop = ''
        if rot[i, 0] == 1:
            symop += 'x'
        elif rot[i, 0] == -1:
            symop += '-x'
        elif rot[i, 0] == 0:
            pass
        else:
            print('error in rotation matrix')
        if rot[i, 1] == 1:
            if symop:
                symop += '+y'
            else:
                symop += 'y'
        elif rot[i, 1] == -1:
            symop += '-y'
        elif rot[i, 1] == 0:
            pass
        else:
            print('error in rotation matrix')
        if rot[i, 2] == 1:
            if symop:
                symop += '+z'
            else:
                symop += 'z'
        elif rot[i, 2] == -1:
            symop += '-z'
        elif rot[i, 2] == 0:
            pass
        else:
            print('error in rotation matrix')
        (shift, dummy) = math.modf(trans[i])
        if shift < 0:
            shift += 1.0
        if abs(shift) > 1.0e-4 and abs(shift - 1.0) > 1.0e-4:
            if shift > 0:
                symop += '+'
            symop += float2rational(trans[i])
        symbols.append(symop)
    return ','.join(symbols)


def roundint(x):
    return int(math.floor(x + 0.5))


def float2rational(x):
    '''make rational representation of x
       assume 0 < x < 1
       Denominator is assumed to be 2,3, or 4, 6.
    '''
    tol = 1.0e-3
    denoms = [2, 3, 4, 6]
    for d in denoms:
        w = x * d
        n = roundint(w)
        if abs(w - n) < tol:
            s = '%d/%d' % (n, d)
            return s


def celllength(cell):
    ''' calculate lattice constants from primitive vectors
        in
           cell[0], cell[1], cell[2]: primitive vectors
        return a,b,c
    '''
    a = cell[0]
    b = cell[1]
    c = cell[2]
    la = math.sqrt(numpy.dot(a, a))
    lb = math.sqrt(numpy.dot(b, b))
    lc = math.sqrt(numpy.dot(c, c))
    return la, lb, lc


def cellangle(cell):
    ''' calculate lattice angles (in degree) from primitive vectors
        in
           cell[0], cell[1], cell[2]: primitive vectors
        return alpha,beta,gamma
    '''
    a = cell[0]
    b = cell[1]
    c = cell[2]
    la, lb, lc = celllength(cell)
    aa = math.degrees(math.acos(numpy.dot(b, c) / (lb * lc)))
    ab = math.degrees(math.acos(numpy.dot(c, a) / (lc * la)))
    ac = math.degrees(math.acos(numpy.dot(a, b) / (la * lb)))
    return aa, ab, ac


def summary(cry, cell, hallnum_ccell, hallnum_pcell, cif, tol=1.0e-5):
    f = open('caldata.txt', 'a')

    # cif
    f.write(cry.matid + ' ')
    chemform = cry.chemical_formula_sum
    chemform = chemform.replace(' ', '')
    f.write(chemform + ' ')
    f.write(cry.space_group_crystal_system + ' ')
    f.write(str(cry.symmetry_Int_Tables_number) + ' ')
    f.write(cry.symmetry_space_group_name + ' ')
    la = cry.cell.length_a
    lb = cry.cell.length_b
    lc = cry.cell.length_c
    aa = cry.cell.angle_alpha
    ab = cry.cell.angle_beta
    ac = cry.cell.angle_gamma
    f.write(
        str(la) + ' ' + str(lb) + ' ' + str(lc) + ' ' + str(aa) + ' ' +
        str(ab) + ' ' + str(ac) + ' ')

    #tol = 1.0e-5

    # conv cell
    sgtype = spglib.get_spacegroup_type(hallnum_ccell)
    spgnum = sgtype['number']
    spgname = sgtype['international_short']
    f.write(str(spgnum) + ' ' + spgname + ' ')

    # prim cell
    sgtype = spglib.get_spacegroup_type(hallnum_pcell)
    spgnum = sgtype['number']
    spgname = sgtype['international_short']
    f.write(str(spgnum) + ' ' + spgname + ' ')

    # calc cell
    scell = spglib.standardize_cell(cell)
    dataset = spglib.get_symmetry_dataset(scell, tol)
    spgnum = dataset['number']
    hallnum = dataset['hall_number']
    sgtype = spglib.get_spacegroup_type(hallnum)
    spgname = sgtype['international_short']
    f.write(str(spgnum) + ' ' + spgname + ' ')
    lattice = scell[0]
    a, b, c = celllength(lattice)
    alpha, beta, gamma = cellangle(lattice)
    f.write(
        str(a) + ' ' + str(b) + ' ' + str(c) + ' ' + str(alpha) + ' ' +
        str(beta) + ' ' + str(gamma) + ' ')

    # cif2cell
    cmd = './cif2cell.sh ' + cif
    os.system(cmd)
    cell, atoms = read_structure_vasp('POSCAR.cif2cell')
    scell = spglib.standardize_cell(cell)
    dataset = spglib.get_symmetry_dataset(scell, tol)
    spgnum = dataset['number']
    hallnum = dataset['hall_number']
    sgtype = spglib.get_spacegroup_type(hallnum)
    spgname = sgtype['international_short']
    f.write(str(spgnum) + ' ' + spgname + ' ')

    f.write('\n')


#-------------------------------------------------------------
if __name__ == '__main__':
    main()
