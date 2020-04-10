#!/bin/env python3

import gemmi

ENANTIOMORPHIC = [76, 78, 91, 92, 95, 96, 144, 145, 151, 152, 153, 154,
                  169, 170, 171, 172, 178, 179, 180, 181, 212, 213]

def generate(spacegroup):
    eu_list = []
    pg_ops = spacegroup.operations()
    pg_ops_list = []
    centering = pg_ops.find_centering()
    for op in pg_ops.sym_ops:
        op.tran = [0, 0, 0]
        pg_ops_list.append(op)
    for sg in gemmi.spacegroup_table():
        ops = sg.operations()
        if ops.find_centering() != centering:
            continue
        ops_list = []
        for op in ops.sym_ops:
            op.tran = [0, 0, 0]
            ops_list.append(op)
        if set(pg_ops_list) == set(ops_list):
            if sg.number in ENANTIOMORPHIC:
                enant_ops = sg.operations()
                enant_ops.change_basis(gemmi.Op('-x,-y,-z'))
                if gemmi.find_spacegroup_by_ops(enant_ops) in eu_list:
                    continue
            eu_list.append(sg)
        if sg.number == 230:
            break
#    print(eu_list)
    for e in eu_list:
        print(e, e.is_reference_setting(), e.number)
    return eu_list

#generate(gemmi.SpaceGroup('P 2 21 2'))
