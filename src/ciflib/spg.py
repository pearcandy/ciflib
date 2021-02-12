import numpy
import spglib


#===================================================
#               standardize cell
#===================================================
def standardize_cell(cry, tol=1.0e-5):
    primcell = numpy.array(cry.unit_prim, float)
    primatom = []

    for a in cry.vatom_prim:
        primatom.append([a.position[0], a.position[1], a.position[2]])
    primnum = []

    for i in range(len(cry.velement_size)):
        for j in range(cry.velement_size[i]):
            primnum.append(i)
    pcell = (primcell, primatom, primnum)
    print_cell(pcell)

    scell = spglib.standardize_cell(pcell)
    dataset = spglib.get_symmetry_dataset(scell, tol)
    symmetry = spglib.get_symmetry(scell)
    print_cell(scell)

    return scell, dataset, symmetry


#===================================================
#                  cif check
#===================================================
def check_by_ciflib(cry, tol=1.0e-5):
    spgnum_cif = cry.symmetry_Int_Tables_number
    convcell = numpy.array(cry.unit_conv, float)
    convatom = []
    for a in cry.vatom_conv:
        convatom.append([a.position[0], a.position[1], a.position[2]])
    convnum = []
    for i in range(len(cry.velement_size_conv)):
        for j in range(cry.velement_size_conv[i]):
            convnum.append(i)

    ccell = (convcell, convatom, convnum)
    #print_cell(ccell)

    dataset = spglib.get_symmetry_dataset(ccell, tol)
    spgnum_ccell = dataset['number']
    hallnum_ccell = dataset['hall_number']
    #print_dataset(dataset)

    cry.symmetry_space_group_name_hall = dataset['hall']

    print('\n[primitive cell]')

    primcell = numpy.array(cry.unit_prim, float)
    primatom = []
    for a in cry.vatom_prim:
        primatom.append([a.position[0], a.position[1], a.position[2]])
    primnum = []
    for i in range(len(cry.velement_size)):
        for j in range(cry.velement_size[i]):
            primnum.append(i)

    pcell = (primcell, primatom, primnum)
    print_cell_name(cry, pcell)

    dataset = spglib.get_symmetry_dataset(pcell, tol)
    spgnum_pcell = dataset['number']
    hallnum_pcell = dataset['hall_number']
    #print_dataset(dataset)

    print('\n[conventional cell]')

    scell = spglib.standardize_cell(pcell)
    print_cell_name(cry, scell)
    dataset = spglib.get_symmetry_dataset(scell, tol)
    spgnum_ccell_pcell_std = dataset['number']
    hallnum_ccell_pcell_std = dataset['hall_number']

    print('\n- check crystal structure')
    print('%-25s %10s %10s' % ('', 'spg num', 'hall num'))
    print('%-25s %10d %10s' % ('cell orginal cif', spgnum_cif, '-'))
    print('%-25s %10d %10d' % ('ciflib> ccell ', spgnum_ccell, hallnum_ccell))
    print('%-25s %10d %10d' % ('ciflib> pcell', spgnum_pcell, hallnum_pcell))
    print('%-25s %10d %10d' %
          ('ccell pcell+std', spgnum_ccell_pcell_std, hallnum_ccell_pcell_std))

    b = True
    if spgnum_cif != spgnum_ccell:
        b = False
        print('wraning: conventional cell... NG')
    if spgnum_cif != spgnum_pcell:
        b = False
        print('wraning: primitive cell... NG')
    if b == True:
        print('ciflib> successfully finished')

    return b, hallnum_ccell, hallnum_pcell


#===================================================
#                  print cell with name index
#===================================================
def print_cell_name(cry, cell):
    lattice, positions, numbers = cell
    print('- ' + 'lattice vectors')
    for vec in lattice:
        print('%15.9f %15.9f %15.9f' % (vec[0], vec[1], vec[2]))
    print('- ' + 'atomic postions and type (frac)')
    for coor, atype in zip(positions, numbers):
        print('%4s %15.9f %15.9f %15.9f' %
              (cry.velement_name[int(atype)], coor[0], coor[1], coor[2]))
    coor_xyz = []
    for coor2 in positions:
        coor_xyz.append(lattice[0] * coor2[0] + lattice[1] * coor2[1] +
                        lattice[2] * coor2[2])
    print('- ' + 'atomic postions and type (xyz)')
    i0 = -1
    for coor, atype in zip(positions, numbers):
        i0 = i0 + 1
        print('%4s %15.9f %15.9f %15.9f' %
              (cry.velement_name[int(atype)], coor_xyz[i0][0], coor_xyz[i0][1],
               coor_xyz[i0][2]))


#===================================================
#                  print cell
#===================================================
def print_cell(cell):
    lattice, positions, numbers = cell
    print('- ' + 'lattice vectors')
    for vec in lattice:
        print('%15.9f %15.9f %15.9f' % (vec[0], vec[1], vec[2]))
    print('- ' + 'atomic postions and type')
    for coor, atype in zip(positions, numbers):
        print('%4d %15.9f %15.9f %15.9f' % (atype, coor[0], coor[1], coor[2]))


#===================================================
#                print dataset
#===================================================
def print_dataset(dataset):
    print('- ' + 'spg_dataset')
    print('  number', dataset['number'])
    print('  international', dataset['international'])
    print('  hall_symbol', dataset['hall'])
    print('  hall_number', dataset['hall_number'])
    print('  choice', dataset['choice'])
    print('  pointgroup', dataset['pointgroup'])
    print('  trans_matrix\n', dataset['transformation_matrix'])
    print('  original_shift\n', dataset['origin_shift'])
    hall_number = dataset['hall_number']
    sgtype = spglib.get_spacegroup_type(hall_number)
    #print_spgtype(sgtype)


#===================================================
#                print spgtype
#===================================================
def print_spgtype(sgtype):
    print('- ' + 'spg_spgtype')
    print('  number', sgtype['number'])
    print('  international_short', sgtype['international_short'])
    print('  international_full', sgtype['international_full'])
    print('  international', sgtype['international'])
    print('  hall_sysmbol', sgtype['hall_symbol'])
    print('  choice', sgtype['choice'])
    print('  schoenflies', sgtype['schoenflies'])
    print('  pointgroup', sgtype['pointgroup_schoenflies'])
    print('  pointgroup_international', sgtype['pointgroup_international'])
    print('  arithmetic_crystal_class_number',
          sgtype['arithmetic_crystal_class_number'])
    print('  arithmetic_crystal_class_symbol',
          sgtype['arithmetic_crystal_class_symbol'])
