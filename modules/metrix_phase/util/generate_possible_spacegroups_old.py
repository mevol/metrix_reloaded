#!/bin/env python3
#
# fast_ep ->
#
# Fast experimental phasing in the spirit of fast_dp, starting from nothing
# and using brute force (and educated guesses) to get everything going.
#
# generate_chiral_spacegroups - generate lists of chiral spacegroups
# for a given pointgroup.
#
# generate_chiral_spacegroups_unique - generate lists of chiral spacegroups
# for a given pointgroup, counting enantiomorphs together (i.e. useful for
# SAD substructure determination.)

#from cctbx.sgtbx import space_group, space_group_symbols, \
#     space_group_symbol_iterator

import gemmi
import re

#def spacegroup_enantiomorph(spacegroup_name):
#    sg = space_group(space_group_symbols(spacegroup_name).hall())
#    enantiomorph = sg.change_basis(sg.type().change_of_hand_op())
#    return enantiomorph.type().lookup_symbol().replace(' ', '')

def spacegroup_enantiomorph(spacegroup_name):
  operations = []
  sg = gemmi.find_spacegroup_by_name(spacegroup_name)
  print(sg)
  sg_ops = sg.operations()
  for op in sg_ops: 
    print(op)
    operations.append(op)
    op_inv = op.inverse()
    print(op_inv)
    operations.append(op_inv)
  print(operations)


#  print(len(sg_ops))
#  print(list(sg_ops.sym_ops))
    
    
#    for op in sg_ops:
#      op_inv = op.inverse()
#      print(op_inv)
    #print(len(sg_ops_inv))
    #print(list(sg_ops_inv.sym_ops))
    
    
    
#    sg_str = re.findall(r'[A-Z]+', str(sg))
#    print(sg_str[-1])
    
    
#    for l in sg_str:
#      if l == "S":
#        pass
#      if l == "G":
#        pass
#      else:  
#        print(l)
    #print(sg_str)
#    pg = sg.point_group_hm()
#    print(pg)
    
#    pg_name = sg_str[-1] + pg
    
#    print(pg_name)
    
#    sys = sg.crystal_system_str()
#    print(sys)
    #for sg in gemmi.spacegroup_table():
      
#    sg_name = gemmi.find_spacegroup_by_name(str(pg_name))
#    print(sg_name)
    
    
#    sg_name_ops = sg_name.operations()
#    print(len(sg_name_ops))
#    print(list(sg_name_ops.sym_ops))
    
#    pg = sg.point_group_hm()
#    print(pg)
#    pg_ops = pg.operations()
#    print(len(pg_ops))
#    print(list(pg_ops.sym_ops))
    #print(list(ops.cen_ops))
#    enantiomorph = gemmi.Op(sg)
#    enantiomorph = _.inverse()
#    print(enantiomorph)
#    return enantiomorph


#def spacegroup_full(spacegroup_name):
#    return space_group(space_group_symbols(
#        spacegroup_name).hall()).type().lookup_symbol()

#def sanitize_spacegroup(spacegroup_name):
#    if not ':' in spacegroup_name:
#        return spacegroup_name
#    assert(spacegroup_name.startswith('R'))
#    return 'H%s' % spacegroup_name.split(':')[0][1:]

#def generate_chiral_spacegroups_unique(pointgroup):
#    sg = space_group(space_group_symbols(pointgroup).hall())
#    pg = sg.build_derived_patterson_group()
#
#    eu_list = []
#
#    for j in space_group_symbol_iterator():
#        sg_test = space_group(j)
#
#        if not sg_test.is_chiral():
#            continue
#
#        pg_test = sg_test.build_derived_patterson_group()
#        if pg_test == pg:
#            enantiomorph = sg_test.change_basis(
#                sg_test.type().change_of_hand_op())
#            if not sg_test in eu_list and not \
#               enantiomorph in eu_list:
#                eu_list.append(sg_test)
#
#    return [
#        sg_test.type().lookup_symbol().replace(' ', '') \
#        for sg_test in eu_list]

#def generate_chiral_spacegroups(pointgroup):
#    sg = space_group(space_group_symbols(pointgroup).hall())
#    pg = sg.build_derived_patterson_group()
#
#    sg_list = []
#
#    for j in space_group_symbol_iterator():
#        sg_test = space_group(j)
#
#        if not sg_test.is_chiral():
#            continue
#
#        pg_test = sg_test.build_derived_patterson_group()
#        if pg_test == pg:
#            if not sg_test in sg_list:
#                sg_list.append(sg_test)
#
#    return [
#        sg_test.type().lookup_symbol().replace(' ', '') \
#        for sg_test in sg_list]

#def test():
#    assert(len(generate_chiral_spacegroups('P422')) == 8)
#    assert(len(generate_chiral_spacegroups_unique('P422')) == 6)
#    assert(len(generate_chiral_spacegroups('P222')) == 8)
#    assert(generate_chiral_spacegroups('P222') == \
#           generate_chiral_spacegroups_unique('P222'))

  print('OK')

if __name__ == '__main__':
    test()
