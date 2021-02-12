#!/usr/bin/env python
import sys
import os
import math
import copy
import collections
import argparse
import re

from . import cell
from . import spgtable
from . import bravais 
from . import symmetry 
from . import atom as atm
from . import position as atom_position
from . import brillouin

def get_float_number(s):
    n = s.find("(")
    #if n > 0:
    #    return float(s[:n])
    #else:
    #    return float(s)
    if n > 0:
        d = float(s[:n])
    else:
        d = float(s)
    sign = 1
    if d < 0:
        d = -d
        sign = -1
    tol = 1.0e-4
    if abs(d-1.0/3.0)<tol:
        d = 1.0/3.0
    if abs(d-2.0/3.0)<tol:
        d = 2.0/3.0
    if abs(d-1.0/6.0)<tol:
        d = 1.0/6.0
    if abs(d-5.0/6.0)<tol:
        d = 5.0/6.0
    return sign * d
        
def parse_line(line):
    blank = " \t"
    squote = "'"
    dquote = '"'
    value_list = []
    word = ""
    in_squote = False
    in_dquote = False
    in_word = False
    length = len(line)
    for i in range(length):
        s = line[i]
        if in_squote:
            if s == squote and (i == length-1 or line[i+1] in blank):
                value_list.append(word)
                word = ""
                in_squote = False
                continue
            else:
                word += s
        elif in_dquote:
            if s == dquote and (i == length-1 or line[i+1] in blank):
                value_list.append(word)
                word = ""
                in_dquote = False
                continue
            else:
                word += s
        elif in_word:
            if s in blank:
                value_list.append(word)
                word = ""
                in_word = False
                continue
            else:
                word += s
        elif s == squote:
            in_squote = True
            word = ""
        elif s == dquote:
            in_dquote = True
            word = ""
        elif s not in blank:
            in_word = True
            word = s
    if len(word) > 0:
        value_list.append(word)
    return value_list


class Crystal():
    version = "0.1.1"

    def __init__( self ):
        self.matid = ""
        self.symmetry_Int_Tables_number = 0
        self.symmetry_cell_setting = ""
        self.space_group_crystal_system = ""
        self.symmetry_space_group_name = ""
        self.symmetry_space_group_name_hall = ""
        self.chemical_formula_sum = ""
        self.publ_section_title = ""
        self.journal_coden_ASTM = ""
        self.journal_name_full = ""
        self.journal_year = ""
        self.journal_volume = ""
        self.journal_page_first = ""
        self.journal_page_last = ""
        self.journal_language = ""
        self.vauthor_name = []

        self.cell = cell.Cell()
        self.spacegroup = spgtable.SpaceGroup()
        self.voperation = []
        self.vatom_cif = []

        self.bravais = bravais.Bravais()
        self.unit_conv = [ atom_position.Position(), atom_position.Position(), atom_position.Position() ]
        self.unit_prim = [ atom_position.Position(), atom_position.Position(), atom_position.Position() ]
        self.recp_conv = [ atom_position.Position(), atom_position.Position(), atom_position.Position() ]
        self.recp_prim = [ atom_position.Position(), atom_position.Position(), atom_position.Position() ]
        self.baseP2C   = [[0,0,0],[0,0,0],[0,0,0]]
        self.baseC2P   = [[0,0,0],[0,0,0],[0,0,0]]

        self.vatom_conv = []
        self.vatom_prim = []

        self.velement_name = []
        self.velement_size = []
        self.velement_size_conv = []

        self.vplane = []
        self.vfacet = []
        self.vsymmL = []

    def read_cif(self, filename):
   
        f = open( filename, 'r' )
        lines = f.readlines()
        f.close()

        self.matid = filename
        n=self.matid.rfind(".")
        if n>=0: self.matid = self.matid[:n]
        n=filename.rfind("-1-2")
        if n>=0: self.matid = self.matid[:n]
        n=filename.rfind("/")
        if n>=0: self.matid = self.matid[n+1:]

        index = 0
        while index < len(lines):
            line = lines[index]
            line = line.strip()
            #print index,line

            m = re.match('(_chemical_formula_sum)\s+\'([^\']+)',line)
            if m: self.chemical_formula_sum = m.group(2)
            m = re.match('(_chemical_name_common)\s+\'([^\']+)',line)
            if m: self.chemical_formula_sum = m.group(2)

            m = re.match('(_space_group_crystal_system)\s+\'([^\']+)',line)
            if m: self.space_group_crystal_system = m.group(2)
            m = re.match('(_space_group_crystal_system)\s+(\w+)',line)
            if m: self.space_group_crystal_system = m.group(2)
            m = re.match('(_symmetry_cell_setting)\s+\'([^\']+)',line)
            if m: self.symmetry_cell_setting = m.group(2)
            m = re.match('(_symmetry_cell_setting)\s+\'(\w+)',line)
            if m: self.symmetry_cell_setting = m.group(2)

            m = re.match('(_symmetry_space_group_name_H-M)\s+\'([^\']+)',line)
            if m: self.symmetry_space_group_name = "".join(m.group(2).split(" "))
            m = re.match('(_symmetry_space_group_name_H-M)\s+\"([^\"]+)',line)
            if m: self.symmetry_space_group_name = "".join(m.group(2).split(" "))
            m = re.match('(_space_group_name_H-M_alt)\s+\'([^\']+)',line)
            if m: self.symmetry_space_group_name = "".join(m.group(2).split(" "))

            m = re.match('(_symmetry_space_group_name_Hall)\s+\'([^\']+)',line)
            if m: self.symmetry_space_group_name_hall = "".join(m.group(2).split(" "))
            m = re.match('(_symmetry_space_group_name_Hall)\s+\"([^\"]+)',line)
            if m: self.symmetry_space_group_name_hall = "".join(m.group(2).split(" "))

            m = re.match('(_symmetry_Int_Tables_number)\s+(\d+)',line)
            if m: self.symmetry_Int_Tables_number = int(m.group(2))
            m = re.match('(_space_group_IT_number)\s+(\d+)',line)
            if m: self.symmetry_Int_Tables_number = int(m.group(2))

            m = re.match('(_cell_length_a)\s+([\d\.]+)',line)
            if m: self.cell.length_a = float(m.group(2))
            m = re.match('(_cell.length_b)\s+([\d\.]+)',line)
            if m: self.cell.length_b = float(m.group(2))
            m = re.match('(_cell.length_c)\s+([\d\.]+)',line)
            if m: self.cell.length_c = float(m.group(2))
            m = re.match('(_cell.angle_alpha)\s+([\d\.]+)',line)
            if m: self.cell.angle_alpha = float(m.group(2))
            m = re.match('(_cell.angle_beta)\s+([\d\.]+)',line)
            if m: self.cell.angle_beta = float(m.group(2))
            m = re.match('(_cell.angle_gamma)\s+([\d\.]+)',line)
            if m: self.cell.angle_gamma = float(m.group(2))

            m = re.match('(_publ_section_title)\s+\'([^\']+)',line)
            if m: self.publ_section_title = m.group(2)
            m = re.match('(_journal_coden_ASTM)\s+(\w+)',line)
            if m: self.journal_coden_ASTM = m.group(2)
            m = re.match('(_journal_name_full)\s+\'([^\']+)',line)
            if m: self.journal_name_full = m.group(2)
            m = re.match('(_journal_year)\s+(\w+)',line)
            if m: self.journal_year = m.group(2)
            m = re.match('(_journal_volume)\s+(\w+)',line)
            if m: self.journal_volume = m.group(2)
            m = re.match('(_journal_page_first)\s+(\w+)',line)
            if m: self.journal_page_first = m.group(2)
            m = re.match('(_journal_page_last)\s+(\w+)',line)
            if m: self.journal_page_last = m.group(2)
            m = re.match('(_journal_language)\s+(\w+)',line)
            if m: self.journal_language = m.group(2)

            if "loop_" == line:
                loop_symm = False
                loop_atom = False
                loop_publ = False
                column       =  0
                column_xyz   = -1
                column_label = -1
                column_x = -1
                column_y = -1
                column_z = -1
                column_occupancy = -1
                column_multiplicity = -1
                column_wyckoff = -1
                column_uiso = -1
                column_symbol = -1
                column_author = -1

                index += 1
                
                while index < len(lines):
                    line = lines[index]
                    line = line.strip()
                    if len(line)==0 or line[0] != "_":
                        break
                    #print (index,line)

                    if False:
                        pass
                    # symmetry
                    elif "_symmetry_equiv_pos_as_xyz" == line:
                        column_xyz = column
                        loop_symm = True
                    elif "_space_group_symop_operation_xyz" == line:
                        column_xyz = column
                        loop_symm = True
                    # atom
                    elif "_atom_site_label" == line:
                        column_label = column
                        loop_atom = True
                    elif "_atom_site_fract_x" == line:
                        column_x = column
                        loop_atom = True
                    elif "_atom_site_fract_y" == line:
                        column_y = column
                        loop_atom = True
                    elif "_atom_site_fract_z" == line:
                        column_z = column
                        loop_atom = True
                    elif "_atom_site_occupancy" == line:
                        column_occupancy = column
                        loop_atom = True
                    elif "_atom_site_symmetry_multiplicity" == line:
                        column_multiplicity = column
                        loop_atom = True
                    elif "_atom_site_Wyckoff_symbol" == line:
                        column_wyckoff = column
                        loop_atom = True
                    elif "_atom_site_U_iso_or_equiv" == line:
                        column_uiso = column
                        loop_atom = True
                    elif "_atom_site_type_symbol" == line:
                        column_symbol = column
                        loop_atom = True
                    # publ
                    elif "_publ_author_name" ==  line:
                        column_author = column
                        loop_publ = True

                    column = column + 1
                    index += 1

                if loop_symm:
                    while index < len(lines):
                        line = lines[index]
                        line = line.strip()
                        if len(line)>0 and line[0]=="#":
                            index += 1
                            continue
                        if len(line)==0 or line == "loop_" or line[0] == "_":
                            index -= 1
                            break
                        #print index,line

                        # terms = line.split()
                        terms = parse_line(line)
                        if len(terms) == 0: break

                        if column_xyz != -1:
                            xyz  = terms[column_xyz]
                            symm = symmetry.SymmetryOperation( xyz )
                            self.voperation.append( symm )

                        index += 1

                if loop_atom:
                    site_index=0
                    while index < len(lines):
                        line = lines[index]
                        line = line.strip()
                        if len(line)>0 and line[0]=="#":
                            index += 1
                            continue
                        if len(line)==0 or line == "loop_" or line[0] == "_":
                            index -= 1
                            break
                        #print (index,line)

                        new_atom = atm.Atom()
                        terms = line.split()
                        if len(terms) == 0: break

                        if column_x != -1:
                            new_atom.position[0] = get_float_number(terms[column_x])
                        if column_y != -1:
                            new_atom.position[1] = get_float_number(terms[column_y])
                        if column_z != -1:
                            new_atom.position[2] = get_float_number(terms[column_z])
                        if column_label != -1:
                            new_atom.label = terms[column_label]
                        if column_occupancy != -1:
                            new_atom.occupancy = get_float_number(terms[column_occupancy])
                        if column_multiplicity != -1:
                            new_atom.multiplicity = int(terms[column_multiplicity])
                        if column_wyckoff != -1:
                            new_atom.wyckoff = terms[column_wyckoff]
                        if column_uiso != -1:
                            new_atom.uiso = terms[column_uiso]
                        if column_symbol != -1:
                            new_atom.element = terms[column_symbol]

                        site_index+=1
                        new_atom.site_index = site_index
                        self.vatom_cif.append(new_atom)

                        index += 1

                if loop_publ:
                    while index < len(lines):
                        line = lines[index]
                        line = line.strip()
                        if len(line)>0 and line[0]=="#":
                            index += 1
                            continue
                        if len(line)==0 or line == "loop_" or line[0] == "_":
                            index -= 1
                            break
                        #print index,line

                        if column_author != -1:
                            m = re.match('\'([^\']+)',line)
                            if m: self.vauthor_name.append(m.group(1))

                        index += 1

            index += 1

    def read_structure_vasp(self,filename):
        f = open( filename, 'r' )

        l = f.readline() ## comment
        l = f.readline() ## factor
        l0 = l.split()
        unit_factor = float(l0[0])

        l = f.readline() ## unitA
        l0 = l.split()
        self.unit_prim[0] = atom_position.Position(float(l0[0]),float(l0[1]),float(l0[2]))
        self.unit_prim[0] *= unit_factor

        l = f.readline() ## unitB
        l0 = l.split()
        self.unit_prim[1] = atom_position.Position(float(l0[0]),float(l0[1]),float(l0[2]))
        self.unit_prim[1] *= unit_factor
        l = f.readline() ## unitC
        l0 = l.split()
        self.unit_prim[2] = atom_position.Position(float(l0[0]),float(l0[1]),float(l0[2]))
        self.unit_prim[2] *= unit_factor

        l = f.readline() ## elements
        l0 = l.split()
        self.velement_name = []
        for l1 in l0:
            self.velement_name.append(l1)

        l = f.readline() ## natom_ele
        l0 = l.split()
        self.velement_size = []
        total = 0
        for l1 in l0:
            self.velement_size.append(int(l1))
            total += int(l1)

        if len(self.velement_name) != len(self.velement_size):
            print(" Error: element name or size.")
            sys.exit(1)

        l = f.readline() ## Direct

        ## without CIF
        if len(self.vatom_prim) == 0:
            for e in range(len(self.velement_name)):
                for a in range(len(self.velement_size)):
                    l = f.readline() ## atoms
                    l0 = l.split()
                    new_atom = atm.Atom()
                    new_atom.position = atom_position.Position(float(l0[0]),float(l0[1]),float(l0[2]))
                    new_atom.label = "%s%d" % (self.velement_name[e], a+1)
                    new_atom.occupancy = 1.0
                    new_atom.multiplicity = 1
                    new_atom.wyckoff = "a"
                    new_atom.uiso = "1.0"
                    new_atom.element = self.velement_name[e]
                    self.vatom_prim.append(new_atom)
        ## with CIF
        elif len(self.vatom_prim) == total:
            for a in range(len(self.vatom_prim)):
                l = f.readline() ## atoms
                l0 = l.split()
                self.vatom_prim[a].position = atom_position.Position(float(l0[0]),float(l0[1]),float(l0[2]))
            for e in range(len(self.velement_name)):
                self.velement_size[e]=0
                self.velement_size_conv[e]=0
        else:
            print(" Error: atom element mismatches.")
            sys.exit(1)

        f.close()

    def write_cif(self,filename):
        f = open(filename,'w')

        f.write("data_\n")
        f.write("\n")
        if self.chemical_formula_sum != "":
            f.write("_chemical_formula_sum '%s'\n" % self.chemical_formula_sum)
        f.write("_space_group_crystal_system '%s'\n" % self.space_group_crystal_system)
        f.write("_symmetry_space_group_name_H-M '%s'\n" % self.symmetry_space_group_name)
        f.write("_symmetry_Int_Tables_number %d\n" % self.symmetry_Int_Tables_number)
        if self.symmetry_space_group_name_hall != "":
            f.write("_symmetry_space_group_name_Hall '%s'\n" % self.symmetry_space_group_name_hall)
        if self.symmetry_cell_setting != "":
            f.write("_symmetry_cell_setting '%s'\n" % self.symmetry_cell_setting)
        f.write("\n")

        f.write("_cell_length_a %f\n" % atom_position.Position_length(self.unit_conv[0]))
        f.write("_cell_length_b %f\n" % atom_position.Position_length(self.unit_conv[1]))
        f.write("_cell_length_c %f\n" % atom_position.Position_length(self.unit_conv[2]))
        f.write("_cell_angle_alpha %f\n" % atom_position.Position_angle(self.unit_conv[1],self.unit_conv[2]))
        f.write("_cell_angle_beta %f\n" % atom_position.Position_angle(self.unit_conv[2],self.unit_conv[0]))
        f.write("_cell_angle_gamma %f\n" % atom_position.Position_angle(self.unit_conv[0],self.unit_conv[1]))
        f.write("\n")

        f.write("loop_\n")
        f.write("_symmetry_equiv_pos_site_id\n")
        f.write("_symmetry_equiv_pos_as_xyz\n")
        if len(self.voperation)>0:
            for i,op in enumerate(self.voperation):
                f.write("%d %s\n" % (i+1,op.codes))
        else:
            f.write( "%d %s\n" % (1,"x,y,z") )
        f.write("\n")

        f.write("loop_\n")
        f.write("_atom_site_label\n")
        f.write("_atom_site_fract_x\n")
        f.write("_atom_site_fract_y\n")
        f.write("_atom_site_fract_z\n")
        f.write("_atom_site_occupancy\n")
        f.write("_atom_site_symmetry_multiplicity\n")
        f.write("_atom_site_Wyckoff_symbol\n")
        f.write("_atom_site_type_symbol\n")
        for a in self.vatom_cif:
            f.write("%s %f %f %f %f %d %s %s\n" % \
                        (a.label,a.position[0],a.position[1],a.position[2], \
                             a.occupancy,a.multiplicity,a.wyckoff,a.element))
        f.write("\n")
        f.close()

    def write_cif2(self, filename):
        f = open(filename,'w')

        f.write("data_\n")
        f.write("\n")
        if self.chemical_formula_sum != "":
            f.write("_chemical_formula_sum '%s'\n" % self.chemical_formula_sum)

        f.write("_space_group_crystal_system '%s'\n" % self.space_group_crystal_system)
        f.write("_symmetry_space_group_name_H-M '%s'\n" % self.symmetry_space_group_name)
        f.write("_symmetry_Int_Tables_number %d\n" % self.symmetry_Int_Tables_number)
        if self.symmetry_cell_setting != "":
            f.write("_symmetry_cell_setting '%s'\n" % self.symmetry_cell_setting)
        f.write("\n")

        f.write("_cell_length_a %f\n" % atom_position.Position_length(self.unit_conv[0]))
        f.write("_cell_length_b %f\n" % atom_position.Position_length(self.unit_conv[1]))
        f.write("_cell_length_c %f\n" % atom_position.Position_length(self.unit_conv[2]))
        f.write("_cell_angle_alpha %f\n" % atom_position.Position_angle(self.unit_conv[1],self.unit_conv[2]))
        f.write("_cell_angle_beta %f\n" % atom_position.Position_angle(self.unit_conv[2],self.unit_conv[0]))
        f.write("_cell_angle_gamma %f\n" % atom_position.Position_angle(self.unit_conv[0],self.unit_conv[1]))
        f.write("\n")

        f.write("loop_\n")
        f.write("_symmetry_equiv_pos_site_id\n")
        f.write("_symmetry_equiv_pos_as_xyz\n")
        if len(self.voperation)>0:
            for i,op in enumerate(self.voperation):
                f.write( "%d %s\n" % (i+1,op.codes) )
        else:
            f.write( "%d %s\n" % (1,"x,y,z") )
        f.write("\n")

        f.write("loop_\n")
        f.write("_atom_site_label\n")
        f.write("_atom_site_fract_x\n")
        f.write("_atom_site_fract_y\n")
        f.write("_atom_site_fract_z\n")
        f.write("_atom_site_occupancy\n")
        f.write("_atom_site_symmetry_multiplicity\n")
        f.write("_atom_site_Wyckoff_symbol\n")
        f.write("_atom_site_type_symbol\n")
        for a in self.vatom_cif:
            f.write("%s %f %f %f %f %d %s %s\n" % \
                        (a.label,a.position[0],a.position[1],a.position[2], \
                             a.occupancy,a.multiplicity,a.wyckoff,a.element))
        f.write("\n")

        f.write("# Conventional cell\n" \
                    "# %12.6f %12.6f %12.6f\n" \
                    "# %12.6f %12.6f %12.6f\n" \
                    "# %12.6f %12.6f %12.6f\n" % \
                    (self.unit_conv[0][0], self.unit_conv[0][1], self.unit_conv[0][2],
                     self.unit_conv[1][0], self.unit_conv[1][1], self.unit_conv[1][2],
                     self.unit_conv[2][0], self.unit_conv[2][1], self.unit_conv[2][2]))

        f.write("# nele:%d natom:%d\n" % (len(self.velement_size),len(self.vatom_conv)))
        f.write("#")
        for e in self.velement_name:
            f.write(" %s" % e)
        f.write("\n")
        f.write("#")
        for e in self.velement_size:
            f.write(" %d" % e)
        f.write("\n")
        for a in self.vatom_conv:
            f.write("#%4s  %f  %f  %f  %d  %d\n" % \
                        (a.element,a.position[0],a.position[1],a.position[2], \
                             a.site_index,a.operation_index))
        f.write("\n")
        f.write("# Primitive cell\n" \
                 "# %12.6f %12.6f %12.6f\n" \
                 "# %12.6f %12.6f %12.6f\n" \
                 "# %12.6f %12.6f %12.6f\n" % \
                    (self.unit_prim[0][0], self.unit_prim[0][1], self.unit_prim[0][2], \
                         self.unit_prim[1][0], self.unit_prim[1][1], self.unit_prim[1][2], \
                         self.unit_prim[2][0], self.unit_prim[2][1], self.unit_prim[2][2]))

        f.write("# nele:%d natom:%d\n" % (len(self.velement_size), len(self.vatom_prim)))
        f.write("#")
        for e in self.velement_name:
            f.write(" %s" % e)
        f.write("\n")
        f.write("#")
        for e in self.velement_size:
            f.write(" %d" % e)
        f.write("\n")
        for a in self.vatom_prim:
            f.write("#%4s  %f  %f  %f  %d  %d\n" % \
                        (a.element,a.position[0],a.position[1],a.position[2], \
                             a.site_index,a.operation_index))
        f.write("\n")

        f.write("_publ_section_title '%s'\n" % self.publ_section_title)
        f.write("_journal_coden_ASTM %s\n" % self.journal_coden_ASTM)
        f.write("_journal_name_full '%s'\n" % self.journal_name_full)
        f.write("_journal_year %s\n" % self.journal_year) 
        f.write("_journal_volume %s\n" % self.journal_volume)
        f.write("_journal_page_first %s\n" % self.journal_page_first)
        f.write("_journal_page_last %s\n" % self.journal_page_last)
        f.write("_journal_language %s\n" % self.journal_language)
        f.write("\n")
        if len(self.vauthor_name)>0:
            f.write("loop_\n")
            f.write("_publ_author_name\n")
            for a in self.vauthor_name:
                f.write("'%s'\n" % a)

        f.close()

    def create_sql(self,filename):
        f = open(filename,'w')

        matid = self.matid
        sqlmat = "(select material_id from cs_material where material='%s')" % matid

        f.write("update compdb.cs_material " \
                    "set num_atom_conv=%d, " \
                    "num_atom_prim=%d " \
                    "where material_id in " \
                    "(select material_id from %s as tmp);\n" % \
                    (len(self.vatom_conv),len(self.vatom_prim), sqlmat))

        f.write("delete from compdb.cs_conventional_cell_atom where material_id=%s;\n" % sqlmat)
        f.write("delete from compdb.cs_conventional_cell where material_id=%s;\n" % sqlmat)

        f.write("insert into compdb.cs_conventional_cell" \
                 "(material_id,cell_a1,cell_a2,cell_a3,cell_b1,cell_b2,cell_b3,cell_c1,cell_c2,cell_c3) " \
                 "values (%s,'%f','%f','%f','%f','%f','%f','%f','%f','%f');\n" % \
                    (sqlmat, \
                         self.unit_conv[0][0], self.unit_conv[0][1], self.unit_conv[0][2], \
                         self.unit_conv[1][0], self.unit_conv[1][1], self.unit_conv[1][2], \
                         self.unit_conv[2][0], self.unit_conv[2][1], self.unit_conv[2][2]))

        f.write("insert into compdb.cs_conventional_cell_atom" \
                     "(material_id,atom_id,type,x,y,z,site_id,op_id) values \n")
        for i,a in enumerate(self.vatom_conv):
            f.write("(%s,%d,'%s','%f','%f','%f',%d,%d)" % \
                        (sqlmat,(i+1), \
                             a.element,a.position[0],a.position[1],a.position[2], \
                             a.site_index,a.operation_index))
            if i<len(self.vatom_conv)-1:
                f.write(",\n")
            else:
                f.write(";\n")

        f.write("delete from compdb.cs_primitive_cell_atom where material_id=%s;\n" % sqlmat)
        f.write("delete from compdb.cs_primitive_cell where material_id=%s;\n" % sqlmat)

        f.write("insert into compdb.cs_primitive_cell" \
                    "(material_id,cell_a1,cell_a2,cell_a3,cell_b1,cell_b2,cell_b3,cell_c1,cell_c2,cell_c3) " \
                    "values (%s,'%f','%f','%f','%f','%f','%f','%f','%f','%f');\n" % \
                    (sqlmat, \
                         self.unit_prim[0][0], self.unit_prim[0][1], self.unit_prim[0][2], \
                         self.unit_prim[1][0], self.unit_prim[1][1], self.unit_prim[1][2], \
                         self.unit_prim[2][0], self.unit_prim[2][1], self.unit_prim[2][2] ) )
        f.write("insert into compdb.cs_primitive_cell_atom" \
                    "(material_id,atom_id,type,x,y,z,site_id,op_id) values \n")
        for i,a in enumerate(self.vatom_prim):
            f.write("(%s,%d,'%s','%f','%f','%f',%d,%d)" % \
                        (sqlmat,(i+1), \
                             a.element,a.position[0],a.position[1],a.position[2], \
                             a.site_index,a.operation_index))
            if i<len(self.vatom_prim)-1:
                f.write(",\n")
            else:
                f.write(";\n")

        f.close()

    def construct_structure(self):
        self.constructVectors()
        self.constructAtoms()
        self.constructBrillouin()

    def construct_structure_p1( self ):
        self.space_group_crystal_system = "triclinic"
        self.symmetry_space_group_name  = "P1"
        self.symmetry_Int_Tables_number = 1
        self.checkBravais2()
        self.recp_conv = atom_position.Position_reciprocal(self.unit_conv)
        self.recp_prim = atom_position.Position_reciprocal(self.unit_prim)
        self.vatom_cif = self.vatom_prim
        self.constructBrillouin()

    def reconstruct_structure(self):
        self.reconstructVectors()

        # update position of atoms in cif
        for a in self.vatom_prim:
            if a.operation_index != 1:
                continue
            self.vatom_cif[a.site_index-1].position = atom_position.Position_translate(self.baseP2C,a.position)
            atom_position.Position_fold(self.vatom_cif[a.site_index-1].position)
        self.constructAtoms()

    def checkBravais2( self ):
        self.bravais.center = ""
        self.cell.angle_alpha=0.0
        self.cell.angle_beta=0.0
        self.cell.angle_gamma=0.0

        eps=1.0e-5

        if self.bravais.center == "" :
            self.unit_conv[0]=self.unit_prim[0]
            self.unit_conv[1]=self.unit_prim[1]
            self.unit_conv[2]=self.unit_prim[2]

            self.cell.angle_alpha = atom_position.Position_angle( self.unit_conv[1], self.unit_conv[2] )
            self.cell.angle_beta = atom_position.Position_angle( self.unit_conv[2], self.unit_conv[0] )
            self.cell.angle_gamma = atom_position.Position_angle( self.unit_conv[0], self.unit_conv[1] )

            if abs(self.cell.angle_alpha-90.00)<eps and \
               abs(self.cell.angle_beta-90.0)<eps and \
               abs(self.cell.angle_gamma-90.0)<eps :
                ## simple
                self.bravais.center = "simple"
            elif abs(self.cell.angle_alpha-90.00)<eps and \
                 abs(self.cell.angle_beta-90.0)<eps and \
                 abs(self.cell.angle_gamma-120.0)<eps :
                ## simple
                self.bravais.center = "simple"

        if self.bravais.center == "" :
            self.unit_conv[0] = - self.unit_prim[0] + self.unit_prim[1] + self.unit_prim[2]
            self.unit_conv[1] = + self.unit_prim[0] - self.unit_prim[1] + self.unit_prim[2]
            self.unit_conv[2] = + self.unit_prim[0] + self.unit_prim[1] - self.unit_prim[2]

            self.cell.angle_alpha = atom_position.Position_angle( self.unit_conv[1], self.unit_conv[2] )
            self.cell.angle_beta = atom_position.Position_angle( self.unit_conv[2], self.unit_conv[0] )
            self.cell.angle_gamma = atom_position.Position_angle( self.unit_conv[0], self.unit_conv[1] )
            if abs(self.cell.angle_alpha-90.00)<eps and \
               abs(self.cell.angle_beta-90.0)<eps and \
               abs(self.cell.angle_gamma-90.0)<eps :
                ## face center
                self.bravais.center = "face"

        if self.bravais.center == "" :
            self.unit_conv[0] = + self.unit_prim[1] + self.unit_prim[2]
            self.unit_conv[1] = + self.unit_prim[2] + self.unit_prim[0]
            self.unit_conv[2] = + self.unit_prim[0] + self.unit_prim[1]

            self.cell.angle_alpha = atom_position.Position_angle( self.unit_conv[1], self.unit_conv[2] )
            self.cell.angle_beta = atom_position.Position_angle( self.unit_conv[2], self.unit_conv[0] )
            self.cell.angle_gamma = atom_position.Position_angle( self.unit_conv[0], self.unit_conv[1] )

            if abs(self.cell.angle_alpha-90.00)<eps and \
               abs(self.cell.angle_beta-90.0)<eps and \
               abs(self.cell.angle_gamma-90.0)<eps :
                ## body center
                self.bravais.center = "body"

        if self.bravais.center == "" :
            self.unit_conv[0] = + self.unit_prim[0] - self.unit_prim[1]
            self.unit_conv[1] = + self.unit_prim[0] + self.unit_prim[1]
            self.unit_conv[2] = + self.unit_prim[2]

            self.cell.angle_alpha = atom_position.Position_angle( self.unit_conv[1], self.unit_conv[2] )
            self.cell.angle_beta = atom_position.Position_angle( self.unit_conv[2], self.unit_conv[0] )
            self.cell.angle_gamma = atom_position.Position_angle( self.unit_conv[0], self.unit_conv[1] )

            if abs(self.cell.angle_alpha-90.00)<eps and \
               abs(self.cell.angle_gamma-90.0)<eps :
                ## base center
                self.bravais.center = "base"

        if self.bravais.center == "" :
            self.unit_conv[0] = self.unit_prim[1] - self.unit_prim[2]
            self.unit_conv[1] = self.unit_prim[2] - self.unit_prim[0]
            self.unit_conv[2] = self.unit_prim[0] + self.unit_prim[1] + self.unit_prim[2]

            self.cell.angle_alpha = atom_position.Position_angle( self.unit_conv[1], self.unit_conv[2] )
            self.cell.angle_beta = atom_position.Position_angle( self.unit_conv[2], self.unit_conv[0] )
            self.cell.angle_gamma = atom_position.Position_angle( self.unit_conv[0], self.unit_conv[1] )

            if abs(self.cell.angle_alpha-90.00)<eps and \
               abs(self.cell.angle_beta-90.0)<eps and \
               abs(self.cell.angle_gamma-120.0)<eps :
                ## trigonal simple
                self.bravais.shape    = "trigonal"
                self.bravais.center = "simple"

        if self.bravais.center == "" :
            self.unit_conv[0]=self.unit_prim[0]
            self.unit_conv[1]=self.unit_prim[1]
            self.unit_conv[2]=self.unit_prim[2]

            self.cell.angle_alpha = atom_position.Position_angle( self.unit_conv[1], self.unit_conv[2] )
            self.cell.angle_beta = atom_position.Position_angle( self.unit_conv[2], self.unit_conv[0] )
            self.cell.angle_gamma = atom_position.Position_angle( self.unit_conv[0], self.unit_conv[1] )

            self.bravais.center = "simple"

        self.cell.length_a = atom_position.Position_length( self.unit_conv[0] )
        self.cell.length_b = atom_position.Position_length( self.unit_conv[1] )
        self.cell.length_c = atom_position.Position_length( self.unit_conv[2] )

        if abs(self.cell.angle_alpha-90.00)<eps and \
           abs(self.cell.angle_beta-90.0)<eps and \
           abs(self.cell.angle_gamma-90.0)<eps :
            ## cubic, tetragonal, orthorhombic
            if abs(self.cell.length_a-self.cell.length_b)<eps and \
               abs(self.cell.length_b-self.cell.length_c)<eps :
                self.bravais.shape = "cubic"
                self.bravais.subtype = "none"
            elif abs(self.cell.length_a-self.cell.length_b)<eps and \
                 self.bravais.center != "face" :
                self.bravais.shape = "tetragonal"
                if self.bravais.center == "simple" :
                    self.bravais.subtype = "none"
                elif self.bravais.center == "body" :
                    if self.cell.length_a < self.cell.length_c :
                        self.bravais.subtype = "type1"
                    else:
                        self.bravais.subtype = "type2"
            else:
                self.bravais.shape = "orthorhombic"
                if self.bravais.center == "simple" :
                    self.bravais.subtype = "none"
                elif self.bravais.center == "face" :
                    ka2 = 1.0/(self.cell.length_a*self.cell.length_a)
                    kb2 = 1.0/(self.cell.length_b*self.cell.length_b)
                    kc2 = 1.0/(self.cell.length_c*self.cell.length_c)
                    if ka2 < kb2 + kc2 and \
                       kb2 < ka2 + kc2 and \
                       kc2 < ka2 + kb2 :
                        self.bravais.subtype = "type1"
                    elif ka2 > kb2 + kc2 :
                        self.bravais.subtype = "type2"
                    elif kb2 > ka2 + kc2 :
                        self.bravais.subtype = "type3"
                    elif kc2 > ka2 + kb2 :
                        self.bravais.subtype = "type4"
                elif self.bravais.center == "body" :
                    if self.cell.length_a>self.cell.length_c and \
                       self.cell.length_a>self.cell.length_b :
                        self.bravais.subtype = "type1"
                    elif self.cell.length_b>self.cell.length_a and \
                         self.cell.length_b>self.cell.length_c :
                        self.bravais.subtype = "type2"
                    elif self.cell.length_c>self.cell.length_a and \
                         self.cell.length_c>self.cell.length_b :
                        self.bravais.subtype = "type3"
                elif self.bravais.center == "base" :
                    ## suppose C base
                    if self.cell.length_a<self.cell.length_b :
                        self.bravais.subtype = "type1"
                    else:
                        self.bravais.subtype = "type2"

        elif abs(self.cell.angle_alpha-90.00)<eps and \
             abs(self.cell.angle_beta-90.0)<eps and \
             abs(self.cell.angle_gamma-120.0)<eps :
            if self.bravais.shape == "trigonal" :
                if self.cell.length_a*math.sqrt(3)<self.cell.length_c*math.sqrt(2) :
                    self.bravais.subtype = "type1"
                else:
                    self.bravais.subtype = "type2"
            else:
                ## hexagonal
                if abs(self.cell.length_a-self.cell.length_b)<eps :
                    self.bravais.shape = "hexagonal"
                    self.bravais.subtype = "none"
                else:
                    self.bravais.shape = "triclinic"
                    self.bravais.subtype = "none"

        elif abs(self.cell.angle_alpha-90.0)<eps and \
             abs(self.cell.angle_gamma-90.0)<eps :
            self.bravais.shape = "monoclinic"
            if self.bravais.center == "simple" :
                self.bravais.subtype = "none"
            elif self.bravais.center == "base" : ## AB plane is the base
                if self.cell.length_a<self.cell.length_b :
                    self.bravais.subtype = "type1"
                else:
                    self.bravais.subtype = "type2"

        elif abs(self.cell.angle_alpha-self.cell.angle_beta)<eps and \
             abs(self.cell.angle_beta-self.cell.angle_gamma)<eps :
            ## trigonal
            if abs(self.cell.length_a-self.cell.length_b)<eps and \
               abs(self.cell.length_b-self.cell.length_c)<eps :
                self.bravais.shape = "trigonal"
                if self.cell.angle_alpha<90.0 :
                    self.bravais.subtype = "type1"
                else:
                    self.bravais.subtype = "type2"
            else:
                self.bravais.shape = "triclinic"

        else:
            self.bravais.shape = "triclinic"

        self.constructVectors()

    def construct_structure_pcell(self):
        try:
            import spg
            cry = Crystal()

            chemform = ""
            num = [0 for i in range(len(self.velement_name))]
            for i,a in enumerate(self.velement_name):
                num[i] = self.velement_size_conv[i]
            if "C" in self.velement_name:
                chemform += "C"+str(num[self.velement_name.index("C")])+" "
            if "H" in self.velement_name:
                chemform += "H"+str(num[self.velement_name.index("H")])+" "
            for i,a in enumerate(self.velement_name):
                if a != "C" and a != "H":
                    chemform += a+str(num[i])+" "
            chemform = chemform.strip()
            cry.chemical_formula_sum = chemform

            scell,dataset,symmetry = spg.standardize_cell(self)            
            hall_number = dataset["hall_number"]
            hall_symbol = dataset["hall"]
            print("hall",hall_number)
            wyckoffs = dataset["wyckoffs"]
            equiv_atoms = dataset["equivalent_atoms"]
            spgnum = dataset["number"]

            cry.symmetry_Int_Tables_number = spgnum
            cry.symmetry_space_group_name = dataset["international"]
            cry.symmetry_space_group_name_hall = dataset["hall"]

            if spgnum < 3:
                crysys = "triclinic"
            elif spgnum < 16:
                crysys =  "monoclinic"
            elif spgnum < 75:
                crysys =  "orthorhombic"
            elif spgnum < 143:
                crysys =  "tetragonal"
            elif spgnum < 168:
                crysys =  "trigonal"
            elif spgnum < 195:
                crysys = "hexagonal"
            elif spgnum < 231:
                crysys = "cubic"
            cry.space_group_crystal_system = crysys

            cry.unit_conv = scell[0]

            for rot,trans in zip(dataset['rotations'], dataset['translations']):
                symbols = make_symop(rot, trans)
                cry.voperation.append(symmetry.SymmetryOperation(symbols))
            
            inequiv_atoms = collections.OrderedDict()
            for a in equiv_atoms:
                if a not in inequiv_atoms:
                    inequiv_atoms[a] = 1
                else:
                    inequiv_atoms[a] += 1

            suffix_list = {}
            for e in self.velement_name:
                suffix_list[e] = 1
            print("ele",self.velement_name)
            print(suffix_list)

            for n in inequiv_atoms:
                new_atom = atm.Atom()
                new_atom.position[0] = self.vatom_prim[n].position[0]
                new_atom.position[1] = self.vatom_prim[n].position[1]
                new_atom.position[2] = self.vatom_prim[n].position[2]
                atom_type = self.vatom_prim[n].element
                new_atom.element = atom_type
                new_atom.occupancy = 1.0
                new_atom.multiplicity = inequiv_atoms[n]
                new_atom.wyckoff = wyckoffs[n]
                new_atom.label = atom_type + str(suffix_list[atom_type])
                suffix_list[atom_type] += 1

                cry.vatom_cif.append(new_atom)

            return cry
        except:
            #import traceback
            #traceback.print_exc()
            return self    

    def write_structure_vasp(self, filename):
        f = open(filename, 'w')

        f.write("%s  %s\n" % (self.matid, self.chemical_formula_sum) )
        f.write("%12.6f\n" % 1.0 )
        f.write("%15.10f%15.10f%15.10f\n" % \
            (self.unit_prim[0][0],self.unit_prim[0][1],self.unit_prim[0][2]) )
        f.write("%15.10f%15.10f%15.10f\n" % \
            (self.unit_prim[1][0],self.unit_prim[1][1],self.unit_prim[1][2]) )
        f.write("%15.10f%15.10f%15.10f\n" % \
            (self.unit_prim[2][0],self.unit_prim[2][1],self.unit_prim[2][2]) )
    
        for e in self.velement_name:
            f.write("%5s" % e)
        f.write("\n")
        for e in self.velement_size:
            f.write("%6d" % e)
        f.write("\n")
        f.write("Direct ! internal coordinates\n")
        for a in self.vatom_prim:
            f.write("%15.10f%15.10f%15.10f  ! %s %d\n" % \
                        (a.position[0],a.position[1],a.position[2],a.element,a.site_index))

        f.close()

    def write_structure_vasp_convcell(self, filename):
        f = open( filename, 'w' )

        f.write("%s  %s\n" % (self.matid, self.chemical_formula_sum) )
        f.write("%12.6f\n" % 1.0 )
        f.write("%12.6f%12.6f%12.6f\n" % \
            (self.unit_conv[0][0],self.unit_conv[0][1],self.unit_conv[0][2]) )
        f.write("%12.6f%12.6f%12.6f\n" % \
            (self.unit_conv[1][0],self.unit_conv[1][1],self.unit_conv[1][2]) )
        f.write("%12.6f%12.6f%12.6f\n" % \
            (self.unit_conv[2][0],self.unit_conv[2][1],self.unit_conv[2][2]) )

        for e in self.velement_name:
            f.write("%5s" % e)
        f.write("\n")
        for e in self.velement_size_conv:
            f.write("%6d" % e)
        f.write("\n")

        f.write("Direct ! internal coordinates\n")
        for a in self.vatom_conv:
            f.write("%12.6f%12.6f%12.6f  ! %s %d\n" % \
                        (a.position[0],a.position[1],a.position[2],a.element,a.site_index))

        f.close()

    def write_structure_espresso(self, filename):
        f = open(filename,'w')

        cfgdir = str(os.path.dirname(__file__))+"/../../config"
        eles = None
        if os.path.exists(cfgdir+"/elements.aw.txt"):
            eles = {}
            f0 = open(cfgdir+"/elements.aw.txt")
            for l in f0:
                if "#" in l:
                    continue
                l0 = l.split()
                atomno = l0[0]
                e = l0[1]
                w = l0[2]
                eles[e] = [atomno,w]

        f.write("\n")
        f.write("CELL_PARAMETERS angstrom\n")
        f.write("%12.6f%12.6f%12.6f\n" % \
            (self.unit_prim[0][0],self.unit_prim[0][1],self.unit_prim[0][2]) )
        f.write("%12.6f%12.6f%12.6f\n" % \
            (self.unit_prim[1][0],self.unit_prim[1][1],self.unit_prim[1][2]) )
        f.write("%12.6f%12.6f%12.6f\n" % \
            (self.unit_prim[2][0],self.unit_prim[2][1],self.unit_prim[2][2]) )
        f.write("\n")        
        f.write("ATOMIC_SPECIES\n")
        if eles != None:
            for e in self.velement_name:
                f.write("%s %s [$PP]\n" % (e,eles[e][1]))
        f.write("ATOMIC_POSITIONS crystal\n")
        for a in self.vatom_prim:
            f.write("%s %12.6f%12.6f%12.6f\n" % (a.element,a.position[0],a.position[1],a.position[2]))
        f.write("\n")
        f.close()

    def write_structure_abinit(self, filename):
        f = open(filename,'w')

        cfgdir = str(os.path.dirname(__file__))+"/../../config"
        eles = None
        if os.path.exists(cfgdir+"/elements.aw.txt"):
            eles = {}
            f0 = open(cfgdir+"/elements.aw.txt")
            for l in f0:
                if "#" in l:
                    continue
                l0 = l.split()
                atomno = l0[0]
                e = l0[1]
                eles[e] = atomno
    
        f.write("\n")
        f.write("acell %12.6f%12.6f%12.6f angstrom\n" % \
                    (self.cell.length_a,self.cell.length_a,self.cell.length_a))
        f.write("rprim\n")
        f.write("%12.6f %12.6f %12.6f\n" % (self.baseP2C[0][0],self.baseP2C[1][0],self.baseP2C[2][0]))
        f.write("%12.6f %12.6f %12.6f\n" % (self.baseP2C[0][1],self.baseP2C[1][1],self.baseP2C[2][1]))
        f.write("%12.6f %12.6f %12.6f\n" % (self.baseP2C[0][2],self.baseP2C[1][2],self.baseP2C[2][2]))
        f.write("\n")
        f.write("ntypat %d\n" % len(self.velement_name))
        f.write("znucl "),
        if eles != None:
            for e in self.velement_name:
                f.write(eles[e]+" ")
        f.write("\n")
        f.write("natom %d\n" % len(self.vatom_prim))
        f.write("typeat ")
        for a in self.vatom_prim:
            for i,e in enumerate(self.velement_name):
                if a.element == e:
                    f.write("%d " % (i+1))
        f.write("\n")           
        f.write("xred\n")
        for a in self.vatom_prim:
            f.write("%12.6f%12.6f%12.6f\n" % (a.position[0],a.position[1],a.position[2]))
        f.write("\n")
        f.close()

    def check_cif( self ):
        spacegrouptable = spgtable.SpaceGroupTable()
        self.spacegroup = spacegrouptable.findByName( self.symmetry_space_group_name )

        if self.spacegroup.number==0:
            print(" Warning: symmetry_space_group_name '%s' in CIF is unknown." % \
              self.symmetry_space_group_name )
            self.spacegroup = spacegrouptable.findByIntTableNumber( self.symmetry_Int_Tables_number )
            if self.spacegroup.number==0:
                print(" Error: symmetry_Int_Tables_number '%d' in CIF is unknown." % \
                  self.symmetry_Int_Tables_number )
                sys.exit(1)

        if False:
            pass
        elif self.space_group_crystal_system == "":
            self.space_group_crystal_system = self.spacegroup.shape
            print(" Warning: space_group_crystal_system not found in CIF, substituted by '%s'" % \
              self.spacegroup.shape )
        elif self.space_group_crystal_system == "trigonal" and self.spacegroup.shape == "hexagonal":
            pass
        elif self.space_group_crystal_system != self.spacegroup.shape:
            print(" Error: space_group_crystal_system '%s' in CIF, mismatches with the shape '%s', detected by nameHM." % \
              (self.space_group_crystal_system, self.spacegroup.shape) )
            sys.exit(1)

        if self.symmetry_Int_Tables_number != self.spacegroup.number:
            print(" Error: symmetry_Int_Tables_number '%d' in CIF, missmatches with the number '%d', detected by nameHM." % \
              (self.symmetry_Int_Tables_number, self.spacegroup.number) )
            sys.exit(1)

        if not self.bravais.check( self.spacegroup, self.cell ):
            print(" Error: broken bravais: %s" % self.spacegroup.shape )
            print("  lenght_a: %f  length_b: %f  length_c: %f" % \
              (self.cell.length_a, self.cell.length_b, self.cell.length_c) )
            print("  angle_alpha: %f  angle_beta: %f  angle_gamma: %f" % \
              (self.cell.angle_alpha, self.cell.angle_beta, self.cell.angle_gamma) )
            return False

        return True

    def print_log( self ):        
        print(" - Crystal Structure:")
        print("    matid: %s" % self.matid)
        print("    chemical_formula: %s" % self.chemical_formula_sum)
        print("    crystal_system: %s" % self.space_group_crystal_system)
        print(" - SpaceGroup:")
        print("   nameHM: %s" % self.spacegroup.name )
        print("   number: %d" % self.spacegroup.number )
        print("   shape : %s" % self.spacegroup.shape  )
        print("   center: %s" % self.spacegroup.center )
        print(" - Bravais Lattice:")
        print("    shape: %s" % self.bravais.shape )
        print("    center: %s" % self.bravais.center )
        print("    subtype: %s" % self.bravais.subtype )
        print(" - Cell:")
        print("    length_a: %f" % self.cell.length_a)
        print("    length_b: %f" % self.cell.length_b)
        print("    length_c: %f" % self.cell.length_c)
        print("    angle_alpha: %f" % self.cell.angle_alpha)
        print("    angle_beta: %f" % self.cell.angle_beta)
        print("    angle_gamma: %f" % self.cell.angle_gamma)
        print

        print(" - Conventional unit vectors: [angstrom]")
        print("    A %12.6f %12.6f %12.6f" % \
          (self.unit_conv[0][0], self.unit_conv[0][1], self.unit_conv[0][2]) )
        print("    B %12.6f %12.6f %12.6f" % \
          (self.unit_conv[1][0], self.unit_conv[1][1], self.unit_conv[1][2]) )
        print("    C %12.6f %12.6f %12.6f" % \
          (self.unit_conv[2][0], self.unit_conv[2][1], self.unit_conv[2][2]) )

        print(" - Primitive self.unit vectors: [angstrom]")
        print("    A %12.6f %12.6f %12.6f" % \
          (self.unit_prim[0][0], self.unit_prim[0][1], self.unit_prim[0][2]) )
        print("    B %12.6f %12.6f %12.6f" % \
          (self.unit_prim[1][0], self.unit_prim[1][1], self.unit_prim[1][2]) )
        print("    C %12.6f %12.6f %12.6f" % \
          (self.unit_prim[2][0], self.unit_prim[2][1], self.unit_prim[2][2]) )

        print(" - Conventional reciprocal vectors: [1/angstrom] (2Pi not multiplied)")
        print("    A %12.6f %12.6f %12.6f" % \
          (self.recp_conv[0][0], self.recp_conv[0][1], self.recp_conv[0][2]) )
        print("    B %12.6f %12.6f %12.6f" % \
          (self.recp_conv[1][0], self.recp_conv[1][1], self.recp_conv[1][2]) )
        print("    C %12.6f %12.6f %12.6f" % \
          (self.recp_conv[2][0], self.recp_conv[2][1], self.recp_conv[2][2]) )

        print(" - Primitive reciprocal vectors: [1/angstrom] (2Pi not multiplied)")
        print("    A %12.6f %12.6f %12.6f" % \
          (self.recp_prim[0][0], self.recp_prim[0][1], self.recp_prim[0][2]) )
        print("    B %12.6f %12.6f %12.6f" % \
          (self.recp_prim[1][0], self.recp_prim[1][1], self.recp_prim[1][2]) )
        print("    C %12.6f %12.6f %12.6f" % \
          (self.recp_prim[2][0], self.recp_prim[2][1], self.recp_prim[2][2]) )

        print(" - Conversion matrix for internal coords: from conv to prim in real space")
        print("%12.6f %12.6f %12.6f" % \
          (self.baseC2P[0][0],self.baseC2P[1][0],self.baseC2P[2][0]) )
        print("%12.6f %12.6f %12.6f" % \
          (self.baseC2P[0][1],self.baseC2P[1][1],self.baseC2P[2][1]) )
        print("%12.6f %12.6f %12.6f" % \
          (self.baseC2P[0][2],self.baseC2P[1][2],self.baseC2P[2][2]) )

        print(" - Conversion matrix for internal coords: from prim to conv in real space")
        print("%12.6f %12.6f %12.6f" % \
          (self.baseP2C[0][0],self.baseP2C[1][0],self.baseP2C[2][0]) )
        print("%12.6f %12.6f %12.6f" % \
          (self.baseP2C[0][1],self.baseP2C[1][1],self.baseP2C[2][1]) )
        print("%12.6f %12.6f %12.6f" % \
          (self.baseP2C[0][2],self.baseP2C[1][2],self.baseP2C[2][2]) )

        print(" - Conversion matrix for internal coords: from conv to prim in recp space")
        print("%12.6f %12.6f %12.6f" % \
          (self.baseP2C[0][0],self.baseP2C[0][1],self.baseP2C[0][2]) )
        print("%12.6f %12.6f %12.6f" % \
          (self.baseP2C[1][0],self.baseP2C[1][1],self.baseP2C[1][2]) )
        print("%12.6f %12.6f %12.6f" % \
          (self.baseP2C[2][0],self.baseP2C[2][1],self.baseP2C[2][2]) )

        print(" - Conversion matrix for internal coords: from prim to conv in recp space")
        print("%12.6f %12.6f %12.6f" % \
          (self.baseC2P[0][0],self.baseC2P[0][1],self.baseC2P[0][2]) )
        print("%12.6f %12.6f %12.6f" % \
          (self.baseC2P[1][0],self.baseC2P[1][1],self.baseC2P[1][2]) )
        print("%12.6f %12.6f %12.6f" % \
          (self.baseC2P[2][0],self.baseC2P[2][1],self.baseC2P[2][2]) )

        print(" - Atoms in conventional cell: [internal coordinates]  %d" % \
          len(self.vatom_conv) )
        print("   (element, a, b, c, site-index, operation-index)")
        for a in range(len(self.vatom_conv)):
            print("%4s  %f  %f  %f  %d  %d" % \
              (self.vatom_conv[a].element, \
               self.vatom_conv[a].position[0], \
               self.vatom_conv[a].position[1], \
               self.vatom_conv[a].position[2], \
               self.vatom_conv[a].site_index, \
               self.vatom_conv[a].operation_index) )

        print(" - Atoms in primitive cell: [internal coordinates]  %d" % \
          len(self.vatom_prim) )
        print("   (element, a, b, c, site-index, operation-index)")
        for a in range(len(self.vatom_prim)):
            print("%4s  %f  %f  %f  %d  %d" % \
              (self.vatom_prim[a].element, \
               self.vatom_prim[a].position[0], \
               self.vatom_prim[a].position[1], \
               self.vatom_prim[a].position[2], \
               self.vatom_prim[a].site_index, \
               self.vatom_prim[a].operation_index) )

        return True

    def constructVectors( self ):
        cosa = math.cos(math.radians(self.cell.angle_alpha))
        sina = math.sin(math.radians(self.cell.angle_alpha))
        cosb = math.cos(math.radians(self.cell.angle_beta ))
        sinb = math.sin(math.radians(self.cell.angle_beta ))
        cosc = math.cos(math.radians(self.cell.angle_gamma))
        sinc = math.sin(math.radians(self.cell.angle_gamma))
        eps  = 1.0e-5

        if self.bravais.shape == "cubic":
            # Cubic
            self.unit_conv[0] = atom_position.Position( self.cell.length_a, 0.0, 0.0 )
            self.unit_conv[1] = atom_position.Position( 0.0, self.cell.length_b, 0.0 )
            self.unit_conv[2] = atom_position.Position( 0.0, 0.0, self.cell.length_c )

            if False:
                pass
            elif self.bravais.center == "simple":
                self.unit_prim[0] = self.unit_conv[0]
                self.unit_prim[1] = self.unit_conv[1]
                self.unit_prim[2] = self.unit_conv[2]

                self.baseP2C[0] = [1.0,0.0,0.0]
                self.baseP2C[1] = [0.0,1.0,0.0]
                self.baseP2C[2] = [0.0,0.0,1.0]

                self.baseC2P[0] = [1.0,0.0,0.0]
                self.baseC2P[1] = [0.0,1.0,0.0]
                self.baseC2P[2] = [0.0,0.0,1.0]

            elif self.bravais.center == "face":
                self.unit_prim[0] = (self.unit_conv[1]+self.unit_conv[2])*0.5
                self.unit_prim[1] = (self.unit_conv[2]+self.unit_conv[0])*0.5
                self.unit_prim[2] = (self.unit_conv[0]+self.unit_conv[1])*0.5

                self.baseP2C[0] = [0.0,0.5,0.5]
                self.baseP2C[1] = [0.5,0.0,0.5]
                self.baseP2C[2] = [0.5,0.5,0.0]

                self.baseC2P[0] = [-1.0,1.0,1.0]
                self.baseC2P[1] = [1.0,-1.0,1.0]
                self.baseC2P[2] = [1.0,1.0,-1.0]

            elif self.bravais.center == "body":
                self.unit_prim[0] = (-self.unit_conv[0]+self.unit_conv[1]+self.unit_conv[2])*0.5
                self.unit_prim[1] = (+self.unit_conv[0]-self.unit_conv[1]+self.unit_conv[2])*0.5
                self.unit_prim[2] = (+self.unit_conv[0]+self.unit_conv[1]-self.unit_conv[2])*0.5

                self.baseP2C[0] = [-0.5,0.5,0.5]
                self.baseP2C[1] = [0.5,-0.5,0.5]
                self.baseP2C[2] = [0.5,0.5,-0.5]

                self.baseC2P[0] = [0.0,1.0,1.0]
                self.baseC2P[1] = [1.0,0.0,1.0]
                self.baseC2P[2] = [1.0,1.0,0.0]

            else:
                print(" Error: unknown cell type. %s %s %s" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        elif self.bravais.shape == "hexagonal":
            # Hexagonal
            # defulat bravais.paxis == "C"
            self.unit_conv[0] = atom_position.Position(self.cell.length_a,0.00,0.00)
            self.unit_conv[1] = atom_position.Position(-0.5*self.cell.length_a,math.sqrt(3)*0.5*self.cell.length_a,0.0)
            self.unit_conv[2] = atom_position.Position(0.00,0.00,self.cell.length_c)

            if self.bravais.center == "simple":
                self.unit_prim[0] = self.unit_conv[0]
                self.unit_prim[1] = self.unit_conv[1]
                self.unit_prim[2] = self.unit_conv[2]

                self.baseP2C[0] = [1.0,0.0,0.0]
                self.baseP2C[1] = [0.0,1.0,0.0]
                self.baseP2C[2] = [0.0,0.0,1.0]

                self.baseC2P[0] = [1.0,0.0,0.0]
                self.baseC2P[1] = [0.0,1.0,0.0]
                self.baseC2P[2] = [0.0,0.0,1.0]
            else:
                print(" Error: unknown cell type. %s %s %s" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        elif self.bravais.shape == "trigonal":
            # Trigonal
            if abs(self.cell.angle_gamma-120.0)<eps: # as hexagonal conventional cell
                self.unit_conv[0] = atom_position.Position(self.cell.length_a, 0.00, 0.00 )
                self.unit_conv[1] = atom_position.Position(-0.5*self.cell.length_a,math.sqrt(3)*0.5*self.cell.length_a, 0.0 )
                self.unit_conv[2] = atom_position.Position(0.00, 0.00, self.cell.length_c )

                if self.bravais.center == "simple":
                    self.unit_prim[0] = 2.0/3*self.unit_conv[0] + 1.0/3*self.unit_conv[1] + 1.0/3*self.unit_conv[2]
                    self.unit_prim[1] = -1.0/3*self.unit_conv[0] + 1.0/3*self.unit_conv[1] + 1.0/3*self.unit_conv[2]
                    self.unit_prim[2] = -1.0/3*self.unit_conv[0] - 2.0/3*self.unit_conv[1] + 1.0/3*self.unit_conv[2]

                    self.baseP2C[0] = [+2.0/3,+1.0/3,1.0/3]
                    self.baseP2C[1] = [-1.0/3,+1.0/3,1.0/3]
                    self.baseP2C[2] = [-1.0/3,-2.0/3,1.0/3]

                    self.baseC2P[0] = [1.0,-1.0, 0.0]
                    self.baseC2P[1] = [0.0,+1.0,-1.0]
                    self.baseC2P[2] = [1.0,+1.0,+1.0]

                else:
                    print(" Error: unknown cell type. %s %s %s" % \
                      (self.bravais.shape, self.bravais.center, self.bravais.subtype ) )

            else: # as trigonal conventional cell
                # vec{A}*vec{A} =            a*a + c*c = self.cell.length_a*self.cell.length_a
                # vec{A}*vec{C} = -0.5*a*a + c*c = self.cell.length_a*self.cell.length_a*cosa
                # then 
                #    -a*a + 2*c*c = self.cell.length_a*self.cell.length_a*cosa*2

                a = math.sqrt(2.0/3.0*(1.0-cosa))*self.cell.length_a
                c = math.sqrt(1.0/3.0*(1.0+2.0*cosa))*self.cell.length_a

                self.unit_prim[0] = atom_position.Position(a,0.0,c)
                self.unit_prim[1] = atom_position.Position(-a*0.5,+math.sqrt(3)*a*0.5,c)
                self.unit_prim[2] = atom_position.Position(-a*0.5,-math.sqrt(3)*a*0.5,c)

                if self.bravais.center == "simple":
                    self.unit_conv[0] = self.unit_prim[0] - self.unit_prim[1]
                    self.unit_conv[1] = self.unit_prim[1] - self.unit_prim[2]
                    self.unit_conv[2] = self.unit_prim[0] + self.unit_prim[1] + self.unit_prim[2]

                    self.baseP2C[0] = [ 2.0/3,+1.0/3,1.0/3]
                    self.baseP2C[1] = [-1.0/3,+1.0/3,1.0/3]
                    self.baseP2C[2] = [-1.0/3,-2.0/3,1.0/3]

                    self.baseC2P[0] = [1.0,-1.0, 0.0]
                    self.baseC2P[1] = [0.0, 1.0,-1.0]
                    self.baseC2P[2] = [1.0, 1.0, 1.0]
                else:
                    print(" Error: unknown cell type. %s %s %s" % \
                      (self.bravais.shape, self.bravais.center, self.bravais.subtype ) )
        elif self.bravais.shape == "tetragonal":
            # Tetragonal
            self.unit_conv[0] = atom_position.Position(self.cell.length_a,0.0,0.0)
            self.unit_conv[1] = atom_position.Position(0.0,self.cell.length_b,0.0)
            self.unit_conv[2] = atom_position.Position(0.0,0.0,self.cell.length_c)

            if self.bravais.center == "simple":
                self.unit_prim[0] = self.unit_conv[0]
                self.unit_prim[1] = self.unit_conv[1]
                self.unit_prim[2] = self.unit_conv[2]

                self.baseP2C[0] = [1.0,0.0,0.0]
                self.baseP2C[1] = [0.0,1.0,0.0]
                self.baseP2C[2] = [0.0,0.0,1.0]

                self.baseC2P[0] = [1.0,0.0,0.0]
                self.baseC2P[1] = [0.0,1.0,0.0]
                self.baseC2P[2] = [0.0,0.0,1.0]

            elif self.bravais.center == "body":
                self.unit_prim[0] = (-self.unit_conv[0]+self.unit_conv[1]+self.unit_conv[2])*0.5
                self.unit_prim[1] = (+self.unit_conv[0]-self.unit_conv[1]+self.unit_conv[2])*0.5
                self.unit_prim[2] = (+self.unit_conv[0]+self.unit_conv[1]-self.unit_conv[2])*0.5

                self.baseP2C[0] = [-0.5,0.5,0.5]
                self.baseP2C[1] = [0.5,-0.5,0.5]
                self.baseP2C[2] = [0.5,0.5,-0.5]

                self.baseC2P[0] = [0.0,1.0,1.0]
                self.baseC2P[1] = [1.0,0.0,1.0]
                self.baseC2P[2] = [1.0,1.0,0.0]

            else:
                print(" Error: unknown cell type. %s %s %s" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        elif self.bravais.shape == "orthorhombic":
            # Orthorhombic
            # default self.bravais.paxis == "C"
            self.unit_conv[0] = atom_position.Position(+self.cell.length_a,0.0,0.0)
            self.unit_conv[1] = atom_position.Position(0.0,+self.cell.length_b,0.0)
            self.unit_conv[2] = atom_position.Position(0.0,0.0,+self.cell.length_c)

            if False:
                pass
            elif self.bravais.center == "simple":
                self.unit_prim[0] = self.unit_conv[0]
                self.unit_prim[1] = self.unit_conv[1]
                self.unit_prim[2] = self.unit_conv[2]

                self.baseP2C[0] = [1.0,0.0,0.0]
                self.baseP2C[1] = [0.0,1.0,0.0]
                self.baseP2C[2] = [0.0,0.0,1.0]

                self.baseC2P[0] = [1.0,0.0,0.0]
                self.baseC2P[1] = [0.0,1.0,0.0]
                self.baseC2P[2] = [0.0,0.0,1.0]

            elif self.bravais.center == "face":
                self.unit_prim[0] = (self.unit_conv[1]+self.unit_conv[2])*0.5
                self.unit_prim[1] = (self.unit_conv[2]+self.unit_conv[0])*0.5
                self.unit_prim[2] = (self.unit_conv[0]+self.unit_conv[1])*0.5

                self.baseP2C[0] = [0.0,0.5,0.5]
                self.baseP2C[1] = [0.5,0.0,0.5]
                self.baseP2C[2] = [0.5,0.5,0.0]

                self.baseC2P[0] = [-1.0,1.0,1.0]
                self.baseC2P[1] = [1.0,-1.0,1.0]
                self.baseC2P[2] = [1.0,1.0,-1.0]

            elif self.bravais.center == "body":
                self.unit_prim[0] = (-self.unit_conv[0]+self.unit_conv[1]+self.unit_conv[2])*0.5
                self.unit_prim[1] = (+self.unit_conv[0]-self.unit_conv[1]+self.unit_conv[2])*0.5
                self.unit_prim[2] = (+self.unit_conv[0]+self.unit_conv[1]-self.unit_conv[2])*0.5

                self.baseP2C[0] = [-0.5,0.5,0.5]
                self.baseP2C[1] = [0.5,-0.5,0.5]
                self.baseP2C[2] = [0.5,0.5,-0.5]

                self.baseC2P[0] = [0.0,1.0,1.0]
                self.baseC2P[1] = [1.0,0.0,1.0]
                self.baseC2P[2] = [1.0,1.0,0.0]
    
            elif self.bravais.center == "base":
                self.unit_prim[0] = (self.unit_conv[1]+self.unit_conv[0])*0.5
                self.unit_prim[1] = (self.unit_conv[1]-self.unit_conv[0])*0.5
                self.unit_prim[2] = (self.unit_conv[2])

                self.baseP2C[0] = [ 0.5,0.5,0.0]
                self.baseP2C[1] = [-0.5,0.5,0.0]
                self.baseP2C[2] = [ 0.0,0.0,1.0]

                self.baseC2P[0] = [1.0,-1.0,0.0]
                self.baseC2P[1] = [1.0, 1.0,0.0]
                self.baseC2P[2] = [0.0, 0.0,1.0]
    
            else:
                print(" Error: unknown cell type. %s %s %s" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        elif self.bravais.shape == "monoclinic":
            # Monoclinic
            # default self.bravais.paxis == "C"
            self.unit_conv[0] = atom_position.Position(self.cell.length_a,0.0,0.0)
            self.unit_conv[1] = atom_position.Position(0.0,self.cell.length_b,0.0)
            self.unit_conv[2] = atom_position.Position(self.cell.length_c*cosb,0.0,self.cell.length_c*sinb)

            if self.bravais.center == "simple":
                self.unit_prim[0] = self.unit_conv[0]
                self.unit_prim[1] = self.unit_conv[1]
                self.unit_prim[2] = self.unit_conv[2]

                self.baseP2C[0] = [1.0,0.0,0.0]
                self.baseP2C[1] = [0.0,1.0,0.0]
                self.baseP2C[2] = [0.0,0.0,1.0]

                self.baseC2P[0] = [1.0,0.0,0.0]
                self.baseC2P[1] = [0.0,1.0,0.0]
                self.baseC2P[2] = [0.0,0.0,1.0]

            elif self.bravais.center == "base":
                self.unit_prim[0] = (self.unit_conv[1]+self.unit_conv[0])*0.5
                self.unit_prim[1] = (self.unit_conv[1]-self.unit_conv[0])*0.5
                self.unit_prim[2] = (self.unit_conv[2])

                self.baseP2C[0] = [ 0.5,0.5,0.0]
                self.baseP2C[1] = [-0.5,0.5,0.0]
                self.baseP2C[2] = [ 0.0,0.0,1.0]

                self.baseC2P[0] = [1.0,-1.0,0.0]
                self.baseC2P[1] = [1.0, 1.0,0.0]
                self.baseC2P[2] = [0.0, 0.0,1.0]

            else:
                print(" Error: unknown cell type. %s %s %s" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        elif self.bravais.shape == "triclinic":
            # Triclinic
            self.unit_conv[0][0] = self.cell.length_a # A_x
            self.unit_conv[0][1] = 0.0 # A_y
            self.unit_conv[0][2] = 0.0 # A_z

            self.unit_conv[1][0] = self.cell.length_b*cosc # B_x
            self.unit_conv[1][1] = self.cell.length_b*sinc # B_y
            self.unit_conv[1][2] = 0.0 # B_z

            self.unit_conv[2][0] = self.cell.length_c*cosb # C_x
            self.unit_conv[2][1] = self.cell.length_c*(cosa-cosb*cosc)/sinc # C_y
            self.unit_conv[2][2] = math.sqrt(self.cell.length_c*self.cell.length_c-self.unit_conv[2][0]*self.unit_conv[2][0]-self.unit_conv[2][1]*self.unit_conv[2][1]) # C_z

            if self.bravais.center == "simple":
                self.unit_prim[0] = self.unit_conv[0]
                self.unit_prim[1] = self.unit_conv[1]
                self.unit_prim[2] = self.unit_conv[2]

                self.baseP2C[0] = [1.0,0.0,0.0]
                self.baseP2C[1] = [0.0,1.0,0.0]
                self.baseP2C[2] = [0.0,0.0,1.0]

                self.baseC2P[0] = [1.0,0.0,0.0]
                self.baseC2P[1] = [0.0,1.0,0.0]
                self.baseC2P[2] = [0.0,0.0,1.0]
    
            else:
                print(" Error: unknown cell type. %s %s %s" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        else:
            print(" Error: unknown cell type. %s %s %s" % \
              (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        self.recp_conv = atom_position.Position_reciprocal(self.unit_conv)
        self.recp_prim = atom_position.Position_reciprocal(self.unit_prim)

        if self.bravais.shape == "triclinic":
            kalpha = atom_position.Position_angle(self.recp_conv[1],self.recp_conv[2])
            kbeta  = atom_position.Position_angle(self.recp_conv[2],self.recp_conv[0])
            kgamma = atom_position.Position_angle(self.recp_conv[0],self.recp_conv[1])

            if kalpha > 90.0 and kbeta > 90.0 and kgamma > 90.0:
                self.bravais.subtype = "type1a"

            if kalpha < 90.0 and kbeta < 90.0 and kgamma < 90.0:
                self.bravais.subtype = "type1b"

            if kalpha > 90.0 and kbeta > 90.0 and abs(kgamma-90.0)<1e-8:
                self.bravais.subtype = "type2a"

            if kalpha < 90.0 and kbeta < 90.0 and abs(kgamma-90.0)<1e-8:
                self.bravais.subtype = "type2b"

        return True

    def reconstructVectors( self ):

        if self.bravais.shape == "cubic" :
            ## Cubic
            if False:
                pass
            elif self.bravais.center == "simple" :
                self.unit_conv[0] = self.unit_prim[0]
                self.unit_conv[1] = self.unit_prim[1]
                self.unit_conv[2] = self.unit_prim[2]
            elif self.bravais.center == "face" :
                self.unit_conv[0] = -self.unit_prim[0]+self.unit_prim[1]+self.unit_prim[2]
                self.unit_conv[1] = +self.unit_prim[0]-self.unit_prim[1]+self.unit_prim[2]
                self.unit_conv[2] = +self.unit_prim[0]+self.unit_prim[1]-self.unit_prim[2]
            elif self.bravais.center == "body" :
                self.unit_conv[0] = self.unit_prim[1]+self.unit_prim[2]
                self.unit_conv[1] = self.unit_prim[2]+self.unit_prim[0]
                self.unit_conv[2] = self.unit_prim[0]+self.unit_prim[1]
            else:
                print(" Error: unknown cell type. %s %s %s\n" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        elif self.bravais.shape == "hexagonal" :
            ## Hexagonal
            if self.bravais.center == "simple" :
                self.unit_conv[0] = self.unit_prim[0]
                self.unit_conv[1] = self.unit_prim[1]
                self.unit_conv[2] = self.unit_prim[2]
            else:
                print(" Error: unknown cell type. %s %s %s\n" % \
                 (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        elif self.bravais.shape == "trigonal" :
            ## Trigonal
            if self.bravais.center == "simple" :
                ## as hexagonal conventional unit cell
                self.unit_conv[0] = self.unit_prim[0] - self.unit_prim[1]
                self.unit_conv[1] = self.unit_prim[1] - self.unit_prim[2]
                self.unit_conv[2] = self.unit_prim[0] + self.unit_prim[1] + self.unit_prim[2]
            else:
                print(" Error: unknown cell type. %s %s %s\n" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        elif self.bravais.shape == "tetragonal" :
            ## Tetragonal
            if False:
                pass
            elif self.bravais.center == "simple" :
                self.unit_conv[0] = self.unit_prim[0]
                self.unit_conv[1] = self.unit_prim[1]
                self.unit_conv[2] = self.unit_prim[2]
            elif self.bravais.center == "body" :
                self.unit_conv[0] = self.unit_prim[1]+self.unit_prim[2]
                self.unit_conv[1] = self.unit_prim[2]+self.unit_prim[0]
                self.unit_conv[2] = self.unit_prim[0]+self.unit_prim[1]
            else:
                print(" Error: unknown cell type. %s %s %s\n" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        elif self.bravais.shape == "orthorhombic" :
            ## Orthorhombic
            if False:
                pass
            elif self.bravais.center == "simple" :
                self.unit_conv[0] = self.unit_prim[0]
                self.unit_conv[1] = self.unit_prim[1]
                self.unit_conv[2] = self.unit_prim[2]
            elif self.bravais.center == "face" :
                self.unit_conv[0] = -self.unit_prim[0]+self.unit_prim[1]+self.unit_prim[2]
                self.unit_conv[1] = +self.unit_prim[0]-self.unit_prim[1]+self.unit_prim[2]
                self.unit_conv[2] = +self.unit_prim[0]+self.unit_prim[1]-self.unit_prim[2]
            elif self.bravais.center == "body" :
                self.unit_conv[0] = self.unit_prim[1]+self.unit_prim[2]
                self.unit_conv[1] = self.unit_prim[2]+self.unit_prim[0]
                self.unit_conv[2] = self.unit_prim[0]+self.unit_prim[1]
            elif self.bravais.center == "base" :
                self.unit_conv[0] = self.unit_prim[0]-self.unit_prim[1]
                self.unit_conv[1] = self.unit_prim[0]+self.unit_prim[1]
                self.unit_conv[2] = self.unit_prim[2]
            else:
                print(" Error: unknown cell type. %s %s %s\n" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        elif self.bravais.shape == "monoclinic" :
            ## Monoclinic
            if False:
                pass
            elif self.bravais.center == "simple" :
                self.unit_conv[0] = self.unit_prim[0]
                self.unit_conv[1] = self.unit_prim[1]
                self.unit_conv[2] = self.unit_prim[2]
            elif self.bravais.center == "base" :
                self.unit_conv[0] = self.unit_prim[0]
                self.unit_conv[1] = +self.unit_prim[1]+self.unit_prim[2]
                self.unit_conv[2] = -self.unit_prim[1]+self.unit_prim[2]
            else:
                print(" Error: unknown cell type. %s %s %s\n" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        elif self.bravais.shape == "triclinic" :
            ## Triclinic
            if self.bravais.center == "simple" :
                self.unit_conv[0] = self.unit_prim[0]
                self.unit_conv[1] = self.unit_prim[1]
                self.unit_conv[2] = self.unit_prim[2]
            else:
                print(" Error: unknown cell type. %s %s %s\n" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        else:
            print(" Error: unknown cell type. %s %s %s\n" % \
              (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        recp_conv = atom_position.Position_reciprocal( self.unit_conv )
        recp_prim = atom_position.Position_reciprocal( self.unit_prim )

    def constructAtoms(self):
        # calc position of atoms in conventional cell
        self.vatom_conv = []
        for a in self.vatom_cif:
            for n,op in enumerate(self.voperation):
                atom = copy.deepcopy(a)
                atom.position = op.operate(a.position)
                atom.operation_index = n+1
                if atm.Atom_find(self.vatom_conv,atom):
                    pass
                else:
                    self.vatom_conv.append(atom)

        # count multiplicity
        for a in self.vatom_cif:
            count=0
            for b in self.vatom_conv:
                if a.site_index == b.site_index:
                    count+=1

            # check multiplicity if available
            if a.multiplicity > 0:
                if a.multiplicity != count:
                    print(" Warning: broken multiplicity. atom site index %d" % a.site_index)
                    a.multiplicity = count
            else:
                a.multiplicity = count

        # calc position of atoms in primitive cell
        self.vatom_prim = []
        for a in self.vatom_conv:
            position_prim = atom_position.Position_translate(self.baseC2P,a.position)
            atom_position.Position_fold(position_prim)

            atom = copy.deepcopy(a)
            atom.position = position_prim

            if not atm.Atom_find(self.vatom_prim,atom):
                self.vatom_prim.append(atom)

        # count size of each elemnet
        self.velement_name = []
        self.velement_size = []
        self.velement_size_conv = []

        for a in self.vatom_prim:
            if not a.element in self.velement_name:
                self.velement_name.append(a.element)
                self.velement_size.append(0)
                self.velement_size_conv.append(0)

        for i,e in enumerate(self.velement_name):
            for a in self.vatom_prim:
                if e == a.element:
                    self.velement_size[i]+=1

        for i,e in enumerate(self.velement_name):
            for a in self.vatom_conv:
                if e == a.element:
                    self.velement_size_conv[i]+=1

        vatom_sorting = []
        for e in self.velement_name:
            for a in self.vatom_prim:
                if e == a.element:
                    vatom_sorting.append(a)

        self.vatom_prim = vatom_sorting

        vatom_sorting = []
        for e in self.velement_name:
            for a in self.vatom_conv:
                if e == a.element:
                    vatom_sorting.append(a)

        self.vatom_conv = vatom_sorting

    def constructBrillouin(self):
        # initial trial points, which make a large box, large enough to involve the zone
        #
        #      6------7
        #     /|     /|
        #    4-+----5 |
        #    | |    | |
        #    | 2----+-3
        #    |/     |/
        #    0------1
        #

        position = [ \
          -self.recp_prim[0]-self.recp_prim[1]-self.recp_prim[2], \
          +self.recp_prim[0]-self.recp_prim[1]-self.recp_prim[2], \
          -self.recp_prim[0]+self.recp_prim[1]-self.recp_prim[2], \
          +self.recp_prim[0]+self.recp_prim[1]-self.recp_prim[2], \
          -self.recp_prim[0]-self.recp_prim[1]+self.recp_prim[2], \
          +self.recp_prim[0]-self.recp_prim[1]+self.recp_prim[2], \
          -self.recp_prim[0]+self.recp_prim[1]+self.recp_prim[2], \
          +self.recp_prim[0]+self.recp_prim[1]+self.recp_prim[2] ]

        self.vfacet = []
        self.vfacet.append( brillouin.BrillouinFacet(position[0],position[1],position[3],position[2]) )
        self.vfacet.append( brillouin.BrillouinFacet(position[2],position[3],position[7],position[6]) )
        self.vfacet.append( brillouin.BrillouinFacet(position[6],position[7],position[5],position[4]) )
        self.vfacet.append( brillouin.BrillouinFacet(position[4],position[5],position[1],position[0]) )
        self.vfacet.append( brillouin.BrillouinFacet(position[1],position[5],position[7],position[3]) )
        self.vfacet.append( brillouin.BrillouinFacet(position[2],position[6],position[4],position[0]) )

        # wigner ceitz planes
        self.vplane = []
        for kx in [-1,0,+1]:
            for ky in [-1,0,+1]:
                for kz in [-1,0,+1]:
                    if kx==0 and ky==0 and kz==0:
                        continue

                    normal = self.recp_prim[0]*(kx*0.5) \
                           + self.recp_prim[1]*(ky*0.5) \
                           + self.recp_prim[2]*(kz*0.5)
                    self.vplane.append(brillouin.BrillouinPlane(normal))
                    
        # cut zone
        for p in range(len(self.vplane)):
            vsegmentc = []

            f=0
            while f < len(self.vfacet):
                if not self.vfacet[f].cut( self.vplane[p] ):
                    self.vfacet.pop( f )
                    continue
                if self.vfacet[f].cutted:
                    if not brillouin.BrillouinSegment_find( vsegmentc, self.vfacet[f].segmentc ):
                        vsegmentc.append( self.vfacet[f].segmentc )
                f+=1

            if len(vsegmentc)>=3:
                facet = brillouin.BrillouinFacet()
                facet.vsegment = copy.deepcopy(vsegmentc)
                facet.order()
                self.vfacet.append( facet )

        # calc brillouin zone
        self.setupKPath()

        return True

    def setupKPath( self ):
        self.vsymmL = []

        if self.bravais.shape == "cubic":
            # Cubic
            if self.bravais.center == "simple":
                # Cubic simple
                self.vsymmL.append( brillouin.BrillouinSegment("Delta", \
                    "Gamma",atom_position.Position(0.0,0.0,0.0), \
                    "X",atom_position.Position(0.5,0.0,0.0) )) # Delta line, G to X
                self.vsymmL.append( brillouin.BrillouinSegment("Z", \
                    "X",atom_position.Position(0.5,0.0,0.0), \
                    "M",atom_position.Position(0.5,0.5,0.0) )) # Z line, X to M
                self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                    "M",atom_position.Position(0.5,0.5,0.0), \
                    "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Sigma line, M to G
                self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                    "Gamma",atom_position.Position(0.0,0.0,0.0), \
                    "R",atom_position.Position(0.5,0.5,0.5) )) # Lambda line, G to R
                self.vsymmL.append( brillouin.BrillouinSegment("T", \
                    "R",atom_position.Position(0.5,0.5,0.5), \
                    "M",atom_position.Position(0.5,0.5,0.0) )) # T line, R to M
                self.vsymmL.append( brillouin.BrillouinSegment("Z", \
                    "M",atom_position.Position(0.5,0.5,0.0), \
                    "X",atom_position.Position(0.5,0.0,0.0) )) # Z line again, M to X
                self.vsymmL.append( brillouin.BrillouinSegment("S", \
                    "X",atom_position.Position(0.5,0.0,0.0), \
                    "R",atom_position.Position(0.5,0.5,0.5) )) # S line, X to R

            elif self.bravais.center == "face":
                # Cubic face
                self.vsymmL.append( brillouin.BrillouinSegment("Delta", \
                    "Gamma",atom_position.Position(0.0,0.0,0.0), \
                    "X",atom_position.Position(1.0,0.0,0.0) )) # Delta line, G to X
                self.vsymmL.append( brillouin.BrillouinSegment("Z", \
                    "X",atom_position.Position(1.0,0.0,0.0), \
                    "W",atom_position.Position(1.0,0.5,0.0) )) # Z line, X to W
                self.vsymmL.append( brillouin.BrillouinSegment("Q", \
                    "W",atom_position.Position(1.0,0.5,0.0), \
                    "L",atom_position.Position(0.5,0.5,0.5) )) # Q line, W to L
                self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                    "L",atom_position.Position(0.5,0.5,0.5), \
                    "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Lambda line, L to G
                self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                    "Gamma",atom_position.Position(0.0,0.0,0.0), \
                    "U/K",atom_position.Position(0.75,0.75,0.0 ) )) # Sigma line, G to U/K
                self.vsymmL.append( brillouin.BrillouinSegment("S", \
                    "U/K",atom_position.Position(1.0,0.5,0.5), \
                    "X",atom_position.Position(1.0,0.0,0.0) )) # S line, U/K to

            elif self.bravais.center == "body":
                # Cubic body
                self.vsymmL.append( brillouin.BrillouinSegment("Delta", \
                    "Gamma",atom_position.Position(0.0,0.0,0.0), \
                    "H",atom_position.Position(1.0,0.0,0.0) )) # Delta line, G to H
                self.vsymmL.append( brillouin.BrillouinSegment("G", \
                    "H",atom_position.Position(1.0,0.0,0.0), \
                    "N",atom_position.Position(0.5,0.5,0.0) )) # G line, H to N
                self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                    "N",atom_position.Position(0.5,0.5,0.0), \
                    "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Sigma line, N to G
                self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                    "Gamma",atom_position.Position(0.0,0.0,0.0), \
                    "P",atom_position.Position(0.5,0.5,0.5) )) # Lambda line, G to P
                self.vsymmL.append( brillouin.BrillouinSegment("D", \
                    "P",atom_position.Position(0.5,0.5,0.5), \
                    "N",atom_position.Position(0.5,0.5,0.0) )) # D line, P to N

                self.vsymmL.append( brillouin.BrillouinSegment("F", \
                    "H",atom_position.Position(1.0,0.0,0.0), \
                    "P",atom_position.Position(0.5,0.5,0.5) )) # F line, H to P

            else:
                print(" Error: unknown cell type. %s %s %s" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        elif self.bravais.shape == "hexagonal":
            # Hexagonal
            self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                "Gamma",atom_position.Position(0.0,0.0,0.0), \
                "M",atom_position.Position(0.5,0.0,0.0) )) # Sigma line, G to M
            self.vsymmL.append( brillouin.BrillouinSegment("T'", \
                "M",atom_position.Position(0.5,0.0,0.0), \
                "K",atom_position.Position(1.0/3,1.0/3,0.0) )) # T' line, M to K
            self.vsymmL.append( brillouin.BrillouinSegment("T", \
                "K",atom_position.Position(1.0/3,1.0/3,0.0), \
                "Gamma",atom_position.Position(0.0,0.0,0.0) )) # T line, K to G
            self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                "Gamma",atom_position.Position(0.0,0.0,0.0), \
                "A",atom_position.Position(0.0,0.0,0.5) )) # Lambda line, G to A
            self.vsymmL.append( brillouin.BrillouinSegment("R", \
                "A",atom_position.Position(0.0,0.0,0.5), \
                "L",atom_position.Position(0.5,0.0,0.5) )) # R line, A to L
            self.vsymmL.append( brillouin.BrillouinSegment("S'", \
                "L",atom_position.Position(0.5,0.0,0.5), \
                "H",atom_position.Position(1.0/3,1.0/3,0.5) )) # S' line, L to H
            self.vsymmL.append( brillouin.BrillouinSegment("S", \
                "H",atom_position.Position(1.0/3,1.0/3,0.5), \
                "A",atom_position.Position(0.0,0.0,0.5) )) # S line, H to A

            self.vsymmL.append( brillouin.BrillouinSegment("U", \
                "L",atom_position.Position(0.5,0.0,0.5), \
                "M",atom_position.Position(0.5,0.0,0.0) )) # U line, L to M

            self.vsymmL.append( brillouin.BrillouinSegment("P", \
                "K",atom_position.Position(1.0/3,1.0/3,0.0), \
                "H",atom_position.Position(1.0/3,1.0/3,0.5) )) # P line, K to H

        elif self.bravais.shape == "trigonal":
            if self.bravais.subtype == "type1":
                # Trigonal type1 theta<90
                eta=1.0
                self.vsymmL.append( brillouin.BrillouinSegment("Q", \
                    "F",atom_position.Position(0.0,0.5,1.0), \
                    "n1",atom_position.Position(2*eta,0.5-eta,1.0) )) # Q line, F to n1
                self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                    "n1",atom_position.Position(eta,eta,0.0), \
                    "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Sigma line, n1 to G
                self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                    "Gamma",atom_position.Position(0.0,0.0,0.0), \
                    "Z",atom_position.Position(0.0,0.0,1.5) )) # Lambda line, G to Z
                self.vsymmL.append( brillouin.BrillouinSegment("B", \
                    "Z",atom_position.Position(0.0,0.0,1.5), \
                    "n2",atom_position.Position(-eta,2*eta,1.5) )) # B line, Z to n2
                self.vsymmL.append( brillouin.BrillouinSegment("Y", \
                    "n2",atom_position.Position(0.5-eta,2*eta,0.5), \
                    "L",atom_position.Position( 0.5,0.0,0.5) )) # Y line, n2 to L
                self.vsymmL.append( brillouin.BrillouinSegment("l1", \
                    "L",atom_position.Position( 0.5,0.0,0.5), \
                    "Gamma",atom_position.Position(0.0,0.0,0.0) )) # L to G

            elif self.bravais.subtype == "type2":
                # Trigonal type2 theta>90
                eta=1.0
                self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                    "F",atom_position.Position(0.5,0.5,0.0), \
                    "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Sigma line, F to G
                self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                    "Gamma",atom_position.Position(0.0,0.0,0.0), \
                    "n1",atom_position.Position(0.0,0.0,eta) )) # Lambda line, G to n1
                self.vsymmL.append( brillouin.BrillouinSegment("P", \
                    "n1",atom_position.Position(1.0,0.0,eta), \
                    "Z",atom_position.Position(1.0,0.0,-0.5) )) # P line, n1 to Z
                self.vsymmL.append( brillouin.BrillouinSegment("Y", \
                    "Z",atom_position.Position(0.0,1.0,0.5), \
                    "L",atom_position.Position(0.5,0.0,0.5) )) # Y line, Z to L
                self.vsymmL.append( brillouin.BrillouinSegment("l1", \
                    "L",atom_position.Position( 0.5,0.0,0.5), \
                    "Gamma",atom_position.Position(0.0,0.0,0.0) )) # L to G
            else:
                print(" Error: unknown cell type. %s %s %s" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        elif self.bravais.shape == "tetragonal":
            if self.bravais.center == "simple":
                # Tetragonal simple
                self.vsymmL.append( brillouin.BrillouinSegment("Delta", \
                    "Gamma",atom_position.Position(0.0,0.0,0.0), \
                    "X",atom_position.Position(0.5,0.0,0.0) )) # Delta line, G to X
                self.vsymmL.append( brillouin.BrillouinSegment("Y", \
                    "X",atom_position.Position(0.5,0.0,0.0), \
                    "M",atom_position.Position(0.5,0.5,0.0) )) # Y line, X to M
                self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                    "M",atom_position.Position(0.5,0.5,0.0), \
                    "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Sigma line, M to G
                self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                    "Gamma",atom_position.Position(0.0,0.0,0.0), \
                    "Z",atom_position.Position(0.0,0.0,0.5) )) # Lambda line, G to Z
                self.vsymmL.append( brillouin.BrillouinSegment("U", \
                    "Z",atom_position.Position(0.0,0.0,0.5), \
                    "R",atom_position.Position(0.5,0.0,0.5) )) # U line, Z to R
                self.vsymmL.append( brillouin.BrillouinSegment("T", \
                    "R",atom_position.Position(0.5,0.0,0.5), \
                    "A",atom_position.Position(0.5,0.5,0.5) )) # T line, R to A
                self.vsymmL.append( brillouin.BrillouinSegment("S", \
                    "A",atom_position.Position(0.5,0.5,0.5), \
                    "Z",atom_position.Position(0.0,0.0,0.5) )) # S line, A to Z

                self.vsymmL.append( brillouin.BrillouinSegment("W", \
                    "R",atom_position.Position(0.5,0.0,0.5), \
                    "X",atom_position.Position(0.5,0.0,0.0) )) # W line, R to X

                self.vsymmL.append( brillouin.BrillouinSegment("V", \
                    "M",atom_position.Position(0.5,0.5,0.0), \
                    "A",atom_position.Position(0.5,0.5,0.5) )) # V line, M to A

            elif self.bravais.center == "body":
                if self.bravais.subtype == "type1":
                    # Tetragonal Body type1
                    eta=1.0
                    self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                        "Z",atom_position.Position(0.0,0.0,1.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Lambda line, Z to G
                    self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "n1",atom_position.Position(eta,0.0,0.0) )) # Sigma line, G to n1
                    self.vsymmL.append( brillouin.BrillouinSegment("F", \
                        "n1",atom_position.Position(eta,0.0,1.0), \
                        "Z",atom_position.Position(0.0,0.0,1.0) )) # F line, n1 to Z
                    self.vsymmL.append( brillouin.BrillouinSegment("U", \
                        "Z",atom_position.Position(0.0,0.0,1.0), \
                        "n2",atom_position.Position(eta,eta,1.0) )) # U line, Z to n2
                    self.vsymmL.append( brillouin.BrillouinSegment("Y", \
                        "n2",atom_position.Position(eta,1.0-eta,0.0), \
                        "X",atom_position.Position(0.5,0.5,0.0) )) # Y line, n2 to X
                    self.vsymmL.append( brillouin.BrillouinSegment("Delta", \
                        "X",atom_position.Position(0.5,0.5,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Delta line, X to G

                    self.vsymmL.append( brillouin.BrillouinSegment("W", \
                        "X",atom_position.Position(0.5,0.5,0.0), \
                        "P",atom_position.Position(0.5,0.5,0.5) )) # W line, X to P
                    self.vsymmL.append( brillouin.BrillouinSegment("Q", \
                        "P",atom_position.Position(0.5,0.5,0.5), \
                        "N",atom_position.Position(0.5,0.0,0.5) )) # Q line, P to N

                elif self.bravais.subtype == "type2":
                    # Tetragonal Body type2
                    eta=1.0
                    self.vsymmL.append( brillouin.BrillouinSegment("V", \
                        "Z",atom_position.Position(1.0,0.0,0.0), \
                        "n1",atom_position.Position(1.0,0.0,eta) )) # V line, Z to n1
                    self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                        "n1",atom_position.Position(0.0,0.0,eta), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Lambda line, n1 to G
                    self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "Z",atom_position.Position(1.0,0.0,0.0) )) # Sigma line, G to Z
                    self.vsymmL.append( brillouin.BrillouinSegment("Y", \
                        "Z",atom_position.Position(1.0,0.0,0.0), \
                        "X",atom_position.Position(0.5,0.5,0.0) )) # Y line, Z to X
                    self.vsymmL.append( brillouin.BrillouinSegment("Delta", \
                        "X",atom_position.Position(0.5,0.5,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Delta line, X to G

                    self.vsymmL.append( brillouin.BrillouinSegment("W", \
                        "X",atom_position.Position(0.5,0.5,0.0), \
                        "P",atom_position.Position(0.5,0.5,0.5) )) # W line, X to P
                    self.vsymmL.append( brillouin.BrillouinSegment("Q", \
                        "P",atom_position.Position(0.5,0.5,0.5), \
                        "N",atom_position.Position(0.5,0.0,0.5) )) # Q line, P to N

                else:
                    print(" Error: unknown cell type. %s %s %s" % \
                      (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

            else:
                print(" Error: unknown cell type. %s %s %s" % \
                      (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        elif self.bravais.shape == "orthorhombic":
            if self.bravais.center == "simple":
                # Orthorhombic simple
                self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                    "Gamma",atom_position.Position(0.0,0.0,0.0), \
                    "X",atom_position.Position(0.5,0.0,0.0) )) # Sigma line, G to X
                self.vsymmL.append( brillouin.BrillouinSegment("D", \
                    "X",atom_position.Position(0.5,0.0,0.0), \
                    "S",atom_position.Position(0.5,0.5,0.0) )) # D line, X to S
                self.vsymmL.append( brillouin.BrillouinSegment("C", \
                    "S",atom_position.Position(0.5,0.5,0.0), \
                    "Y",atom_position.Position(0.0,0.5,0.0) )) # C line, S to Y
                self.vsymmL.append( brillouin.BrillouinSegment("Delta", \
                    "Y",atom_position.Position(0.0,0.5,0.0), \
                    "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Delta line, Y to G
                self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                    "Gamma",atom_position.Position(0.0,0.0,0.0), \
                    "Z",atom_position.Position(0.0,0.0,0.5) )) # Lambda line, G to Z
                self.vsymmL.append( brillouin.BrillouinSegment("A", \
                    "Z",atom_position.Position(0.0,0.0,0.5), \
                    "U",atom_position.Position(0.5,0.0,0.5) )) # A line, Z to U
                self.vsymmL.append( brillouin.BrillouinSegment("P", \
                    "U",atom_position.Position(0.5,0.0,0.5), \
                    "R",atom_position.Position(0.5,0.5,0.5) )) # P line, U to R
                self.vsymmL.append( brillouin.BrillouinSegment("E", \
                    "R",atom_position.Position(0.5,0.5,0.5), \
                    "T",atom_position.Position(0.0,0.5,0.5) )) # E line, R to T
                self.vsymmL.append( brillouin.BrillouinSegment("B", \
                    "T",atom_position.Position(0.0,0.5,0.5), \
                    "Z",atom_position.Position(0.0,0.0,0.5) )) # B line, T to Z

                self.vsymmL.append( brillouin.BrillouinSegment("G", \
                    "U",atom_position.Position(0.5,0.0,0.5), \
                    "X",atom_position.Position(0.5,0.0,0.0) )) # G line, U to X

                self.vsymmL.append( brillouin.BrillouinSegment("Q", \
                    "S",atom_position.Position(0.5,0.5,0.0), \
                    "R",atom_position.Position(0.5,0.5,0.5) )) # Q line, S to R

                self.vsymmL.append( brillouin.BrillouinSegment("H", \
                    "T",atom_position.Position(0.0,0.5,0.5), \
                    "Y",atom_position.Position(0.0,0.5,0.0) )) # H line, T to Y

            elif self.bravais.center == "face":
                if self.bravais.subtype == "type1":
                    # Orthorhombic face center type1
                    eta=1.0
                    self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                        "Z",atom_position.Position(0.0,0.0,1.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Lambda line, Z to G
                    self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "X",atom_position.Position(1.0,0.0,0.0) )) # Sigma line, G to X
                    self.vsymmL.append( brillouin.BrillouinSegment("D", \
                        "X",atom_position.Position(1.0,0.0,0.0), \
                        "n1",atom_position.Position(1.0,eta,0.0) )) # D line, X to n1
                    self.vsymmL.append( brillouin.BrillouinSegment("B", \
                        "n1",atom_position.Position(0.0,eta,1.0), \
                        "Z",atom_position.Position(0.0,0.0,1.0) )) # B line, n1 to Z
                    self.vsymmL.append( brillouin.BrillouinSegment("A", \
                        "Z",atom_position.Position(0.0,0.0,1.0), \
                        "n2",atom_position.Position(eta,0.0,1.0) )) # A line, Z to n2
                    self.vsymmL.append( brillouin.BrillouinSegment("C", \
                        "n2",atom_position.Position(eta,1.0,0.0), \
                        "Y",atom_position.Position(0.0,1.0,0.0) )) # C line, n2 to Y
                    self.vsymmL.append( brillouin.BrillouinSegment("Delta", \
                        "Y",atom_position.Position(0.0,1.0,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Delta line, Y to G

                    self.vsymmL.append( brillouin.BrillouinSegment("H", \
                        "Y",atom_position.Position(0.0,1.0,0.0), \
                        "n3",atom_position.Position(0.0,1.0,eta) )) # H line, Y to n3
                    self.vsymmL.append( brillouin.BrillouinSegment("G", \
                        "n3",atom_position.Position(1.0,0.0,eta), \
                        "X",atom_position.Position(1.0,0.0,0.0) )) # G line, n3 to X

                elif self.bravais.subtype == "type2":
                    # Orthorhombic face center type2
                    eta=1.0
                    self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                        "Z",atom_position.Position(0.0,0.0,1.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Lambda line, Z to G
                    self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "n1",atom_position.Position(eta,0.0,0.0) )) # Sigma line, G to n1
                    self.vsymmL.append( brillouin.BrillouinSegment("U", \
                        "n1",atom_position.Position(eta,1.0,1.0), \
                        "X",atom_position.Position(0.0,1.0,1.0) )) # U line, n1 to X
                    self.vsymmL.append( brillouin.BrillouinSegment("B", \
                        "X",atom_position.Position(0.0,1.0,1.0), \
                        "Z",atom_position.Position(0.0,0.0,1.0) )) # B line, X to Z
                    self.vsymmL.append( brillouin.BrillouinSegment("A", \
                        "Z",atom_position.Position(0.0,0.0,1.0), \
                        "n2",atom_position.Position(eta,0.0,1.0) )) # A line, Z to n2
                    self.vsymmL.append( brillouin.BrillouinSegment("C", \
                        "n2",atom_position.Position(eta,1.0,0.0), \
                        "Y",atom_position.Position(0.0,1.0,0.0) )) # C line, n2 to Y
                    self.vsymmL.append( brillouin.BrillouinSegment("Delta", \
                        "Y",atom_position.Position(0.0,1.0,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Delta line, Y to G

                    self.vsymmL.append( brillouin.BrillouinSegment("H", \
                        "Y",atom_position.Position(0.0,1.0,0.0), \
                        "X",atom_position.Position(0.0,1.0,1.0) )) # H line, Y to X

                elif self.bravais.subtype == "type3":
                    # Orthorhombic face center type3
                    eta=1.0
                    self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                        "Z",atom_position.Position(0.0,0.0,1.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Lambda line, Z to G
                    self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "X",atom_position.Position(1.0,0.0,0.0) )) # Sigma line, G to X
                    self.vsymmL.append( brillouin.BrillouinSegment("D", \
                        "X",atom_position.Position(1.0,0.0,0.0), \
                        "n2",atom_position.Position(1.0,eta,0.0) )) # D line, X to n2
                    self.vsymmL.append( brillouin.BrillouinSegment("B", \
                        "n2",atom_position.Position(0.0,eta,1.0), \
                        "Z",atom_position.Position(0.0,0.0,1.0) )) # B line, n2 to Z
                    self.vsymmL.append( brillouin.BrillouinSegment("A", \
                        "Z",atom_position.Position(0.0,0.0,1.0), \
                        "Y",atom_position.Position(1.0,0.0,1.0) )) # A line, Z to Y
                    self.vsymmL.append( brillouin.BrillouinSegment("R", \
                        "Y",atom_position.Position(1.0,0.0,1.0), \
                        "n1",atom_position.Position(1.0,eta,1.0) )) # R line, Y to n1
                    self.vsymmL.append( brillouin.BrillouinSegment("Delta", \
                        "n1",atom_position.Position(0.0,eta,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Delta line, n1 to G

                    self.vsymmL.append( brillouin.BrillouinSegment("G", \
                        "Y",atom_position.Position(1.0,0.0,1.0), \
                        "X",atom_position.Position(1.0,0.0,0.0) )) # G line, Y to X

                elif self.bravais.subtype == "type4":
                    # Orthorhombic face center type4
                    eta=1.0
                    self.vsymmL.append( brillouin.BrillouinSegment("Q", \
                        "Z",atom_position.Position(1.0,1.0,0.0), \
                        "n1",atom_position.Position(1.0,1.0,eta) )) # Q line, Z to n1
                    self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                        "n1",atom_position.Position(0.0,0.0,eta), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Lambda line, n1 to G
                    self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "X",atom_position.Position(1.0,0.0,0.0) )) # Sigma line, G to X
                    self.vsymmL.append( brillouin.BrillouinSegment("D", \
                        "X",atom_position.Position(1.0,0.0,0.0), \
                        "Z",atom_position.Position(1.0,1.0,0.0) )) # D line, X to Z
                    self.vsymmL.append( brillouin.BrillouinSegment("C", \
                        "Z",atom_position.Position(1.0,1.0,0.0), \
                        "Y",atom_position.Position(0.0,1.0,0.0) )) # C line, Z to Y
                    self.vsymmL.append( brillouin.BrillouinSegment("Delta", \
                        "Y",atom_position.Position(0.0,1.0,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Delta line, Y to G

                    self.vsymmL.append( brillouin.BrillouinSegment("H", \
                        "Y",atom_position.Position(0.0,1.0,0.0), \
                        "n2",atom_position.Position(0.0,1.0,eta) )) # H line, Y to n2
                    self.vsymmL.append( brillouin.BrillouinSegment("G", \
                        "n2",atom_position.Position(1.0,0.0,eta), \
                        "X",atom_position.Position(1.0,0.0,0.0) )) # G line, n2 to X

                else:
                    print("Error: unknown cell type. %s %s %s" % \
                      (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

            elif self.bravais.center == "body":
                if self.bravais.subtype == "type1":
                    # Orthorhombic body center type1
                    eta=1.0
                    self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "X",atom_position.Position(1.0,0.0,0.0) )) # Sigma line, G to X
                    self.vsymmL.append( brillouin.BrillouinSegment("U", \
                        "X",atom_position.Position(1.0,0.0,0.0), \
                        "n1",atom_position.Position(1.0,eta,0.0) )) # U line, X to n1
                    self.vsymmL.append( brillouin.BrillouinSegment("Delta", \
                        "n1",atom_position.Position(0.0,eta,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Delta line, n1 to G
                    self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "n2",atom_position.Position(0.0,0.0,eta) )) # Lambda line, G to n2
                    self.vsymmL.append( brillouin.BrillouinSegment("G", \
                        "n2",atom_position.Position(1.0,0.0,eta), \
                        "X",atom_position.Position(1.0,0.0,0.0) )) # G line, n2 to X

                    self.vsymmL.append( brillouin.BrillouinSegment("Q", \
                        "R",atom_position.Position(0.5,0.0,0.5), \
                        "W",atom_position.Position(0.5,0.5,0.5) )) # Q line, R to W
                    self.vsymmL.append( brillouin.BrillouinSegment("D", \
                        "W",atom_position.Position(0.5,0.5,0.5), \
                        "S",atom_position.Position(0.0,0.5,0.5) )) # D line, W to S

                    self.vsymmL.append( brillouin.BrillouinSegment("P", \
                        "T",atom_position.Position(0.5,0.5,0.0), \
                        "W",atom_position.Position(0.5,0.5,0.5) )) # P line, T to W

                elif self.bravais.subtype == "type2":
                    # Orthorhombic body center type2
                    eta=1.0
                    self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "n2",atom_position.Position(eta,0.0,0.0) )) # Sigma line, G to n2
                    self.vsymmL.append( brillouin.BrillouinSegment("F", \
                        "n2",atom_position.Position(eta,1.0,0.0), \
                        "X",atom_position.Position(0.0,1.0,0.0) )) # F line, n2 to X
                    self.vsymmL.append( brillouin.BrillouinSegment("Delta", \
                        "X",atom_position.Position(0.0,1.0,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Delta line, X to G
                    self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "n1",atom_position.Position(0.0,0.0,eta) )) # Lambda line, G to n1
                    self.vsymmL.append( brillouin.BrillouinSegment("G", \
                        "n1",atom_position.Position(0.0,1.0,eta), \
                        "X",atom_position.Position(0.0,1.0,0.0) )) # G line, n1 to X

                    self.vsymmL.append( brillouin.BrillouinSegment("Q", \
                        "R",atom_position.Position(0.5,0.0,0.5), \
                        "W",atom_position.Position(0.5,0.5,0.5) )) # Q line, R to W
                    self.vsymmL.append( brillouin.BrillouinSegment("D", \
                        "W",atom_position.Position(0.5,0.5,0.5), \
                        "S",atom_position.Position(0.0,0.5,0.5) )) # D line, W to S

                    self.vsymmL.append( brillouin.BrillouinSegment("P", \
                        "T",atom_position.Position(0.5,0.5,0.0), \
                        "W",atom_position.Position(0.5,0.5,0.5) )) # P line, T to W

                elif self.bravais.subtype == "type3":
                    # Orthorhombic body center type3
                    eta=1.0
                    self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "n1",atom_position.Position(eta,0.0,0.0) )) # Sigma line, G to n1
                    self.vsymmL.append( brillouin.BrillouinSegment("F", \
                        "n1",atom_position.Position(eta,0.0,1.0), \
                        "X",atom_position.Position(0.0,0.0,1.0) )) # F line, n1 to X
                    self.vsymmL.append( brillouin.BrillouinSegment("U", \
                        "X",atom_position.Position(0.0,0.0,1.0), \
                        "n2",atom_position.Position(0.0,eta,1.0) )) # U line, X to n2
                    self.vsymmL.append( brillouin.BrillouinSegment("Delta", \
                        "n2",atom_position.Position(0.0,eta,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Delta line, n2 to G
                    self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "X",atom_position.Position(0.0,0.0,1.0) )) # Lamnda line, G to X

                    self.vsymmL.append( brillouin.BrillouinSegment("Q", \
                        "R",atom_position.Position(0.5,0.0,0.5), \
                        "W",atom_position.Position(0.5,0.5,0.5) )) # Q line, R to W
                    self.vsymmL.append( brillouin.BrillouinSegment("D", \
                        "W",atom_position.Position(0.5,0.5,0.5), \
                        "S",atom_position.Position(0.0,0.5,0.5) )) # D line, W to S

                    self.vsymmL.append( brillouin.BrillouinSegment("P", \
                        "T",atom_position.Position(0.5,0.5,0.0), \
                        "W",atom_position.Position(0.5,0.5,0.5) )) # P line, T to W

                else:
                    print("Error: unknown cell type. %s %s %s" % \
                      (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

            elif self.bravais.center == "base":
                if self.bravais.subtype == "type1":
                    # Orthorhombic baseC center type1
                    eta=1.0
                    self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "n1",atom_position.Position(eta,0.0,0.0) )) # Sigma line, G to n1
                    self.vsymmL.append( brillouin.BrillouinSegment("C", \
                        "n1",atom_position.Position(eta,1.0,0.0), \
                        "Y",atom_position.Position(0.0,1.0,0.0) )) # C line, n1 to Y
                    self.vsymmL.append( brillouin.BrillouinSegment("Delta", \
                        "Y",atom_position.Position(0.0,1.0,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Delta line, Y to G
                    self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "Z",atom_position.Position(0.0,0.0,0.5) )) # Lambda line, G to Z
                    self.vsymmL.append( brillouin.BrillouinSegment("A", \
                        "Z",atom_position.Position(0.0,0.0,0.5), \
                        "n2",atom_position.Position(eta,0.0,0.5) )) # A line, Z to n1
                    self.vsymmL.append( brillouin.BrillouinSegment("E", \
                        "n2",atom_position.Position(eta,1.0,0.5), \
                        "T",atom_position.Position(0.0,1.0,0.5) )) # E line, n1 to T
                    self.vsymmL.append( brillouin.BrillouinSegment("B", \
                        "T",atom_position.Position(0.0,1.0,0.5), \
                        "Z",atom_position.Position(0.0,0.0,0.5) )) # B line, T to Z

                    self.vsymmL.append( brillouin.BrillouinSegment("H", \
                        "Y",atom_position.Position(0.0,1.0,0.0), \
                        "T",atom_position.Position(0.0,1.0,0.5) )) # H line, Y to T

                    self.vsymmL.append( brillouin.BrillouinSegment("D", \
                        "S",atom_position.Position(0.5,0.5,0.0), \
                        "R",atom_position.Position(0.5,0.5,0.5) ))# D line, S to R

                elif self.bravais.subtype == "type2":
                    # Orthorhombic baseC center type2
                    eta=1.0
                    self.vsymmL.append( brillouin.BrillouinSegment("Sigma", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "Y",atom_position.Position(1.0,0.0,0.0) )) # Sigma line, G to Y
                    self.vsymmL.append( brillouin.BrillouinSegment("F", \
                        "Y",atom_position.Position(1.0,0.0,0.0), \
                        "n1",atom_position.Position(1.0,eta,0.0) )) # F line, Y to n1
                    self.vsymmL.append( brillouin.BrillouinSegment("Delta", \
                        "n1",atom_position.Position(0.0,eta,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Delta line, n1 to G
                    self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "Z",atom_position.Position(0.0,0.0,0.5) )) # Lambda line, G to Z
                    self.vsymmL.append( brillouin.BrillouinSegment("A", \
                        "Z",atom_position.Position(0.0,0.0,0.5), \
                        "T",atom_position.Position(1.0,0.0,0.5) )) # A line, Z to T
                    self.vsymmL.append( brillouin.BrillouinSegment("G", \
                        "T",atom_position.Position(1.0,0.0,0.5), \
                        "n2",atom_position.Position(1.0,eta,0.5) )) # G line, T to n2
                    self.vsymmL.append( brillouin.BrillouinSegment("B", \
                        "n2",atom_position.Position(0.0,eta,0.5), \
                        "Z",atom_position.Position(0.0,0.0,0.5) )) # B line, n2 to Z

                    self.vsymmL.append( brillouin.BrillouinSegment("H", \
                        "Y",atom_position.Position(1.0,0.0,0.0), \
                        "T",atom_position.Position(1.0,0.0,0.5) )) # H line, Y to T

                    self.vsymmL.append( brillouin.BrillouinSegment("D", \
                        "S",atom_position.Position(0.5,0.5,0.0), \
                        "R",atom_position.Position(0.5,0.5,0.5) )) # D line, S to R
            else:
                print("Error: unknown cell type. %s %s %s" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        elif self.bravais.shape == "monoclinic":
            if self.bravais.center == "simple":
                # Monoclinic simple
                self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                    "Gamma",atom_position.Position(0.0,0.0,0.0), \
                    "Y",atom_position.Position(0.0,0.5,0.0) )) # Lambda line, G to Y
                self.vsymmL.append( brillouin.BrillouinSegment("l1", \
                    "Y",atom_position.Position(0.0,0.5,0.0), \
                    "C",atom_position.Position(0.0,0.5,0.5) )) # Y to C
                self.vsymmL.append( brillouin.BrillouinSegment("W", \
                    "C",atom_position.Position(0.0,0.5,0.5), \
                    "Z",atom_position.Position(0.0,0.0,0.5) )) # W line, C to Z
                self.vsymmL.append( brillouin.BrillouinSegment("l2", \
                    "Z",atom_position.Position(0.0,0.0,0.5), \
                    "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Z to G
                self.vsymmL.append( brillouin.BrillouinSegment("l3", \
                    "Gamma",atom_position.Position(0.0,0.0,0.0), \
                    "X",atom_position.Position(0.5,0.0,0.0) )) # G to X
                self.vsymmL.append( brillouin.BrillouinSegment("V", \
                    "X",atom_position.Position(0.5,0.0,0.0), \
                    "A",atom_position.Position(0.5,0.5,0.0) )) # V line X to A
                self.vsymmL.append( brillouin.BrillouinSegment("l4", \
                    "A",atom_position.Position(0.5,0.5,0.0), \
                    "Y",atom_position.Position(0.0,0.5,0.0) )) # A to Y
                self.vsymmL.append( brillouin.BrillouinSegment("l5", \
                    "Y",atom_position.Position(0.0,0.5,0.0), \
                    "E",atom_position.Position(0.5,0.5,-0.5) )) # Y to E
                self.vsymmL.append( brillouin.BrillouinSegment("U", \
                    "E",atom_position.Position(0.5,0.5,-0.5), \
                    "D",atom_position.Position(0.5,0.0,-0.5) )) # U line, E to D
                self.vsymmL.append( brillouin.BrillouinSegment("l6", \
                    "D",atom_position.Position(0.5,0.0,-0.5), \
                    "Gamma",atom_position.Position(0.0,0.0,0.0) )) # D to G

            elif self.bravais.center == "base":
                if self.bravais.subtype == "type1":
                    # Monoclinic base center type1
                    eta=1.0
                    self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "Y",atom_position.Position(0.0,1.0,0.0) )) # Lambda line, G to Y
                    self.vsymmL.append( brillouin.BrillouinSegment("l1", \
                        "Y",atom_position.Position(0.0,1.0,0.0), \
                        "n1",atom_position.Position(eta,1.0,0.0) )) # Y to n1
                    self.vsymmL.append( brillouin.BrillouinSegment("l2", \
                        "n1",atom_position.Position(eta,0.0,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) #  n1 to G
                    self.vsymmL.append( brillouin.BrillouinSegment("l3", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "Z",atom_position.Position(0.0,0.0,0.5) )) #  G to Z
                    self.vsymmL.append( brillouin.BrillouinSegment("U", \
                        "Z",atom_position.Position(0.0,0.0,0.5), \
                        "M",atom_position.Position(0.0,1.0,0.5) )) # U line, Z to M
                    self.vsymmL.append( brillouin.BrillouinSegment("l4", \
                        "M",atom_position.Position(0.0,1.0,0.5), \
                        "Y",atom_position.Position(0.0,1.0,0.0) )) # M to Y

                elif self.bravais.subtype == "type2":
                    # Monoclinic baseC center type2
                    eta=1.0
                    self.vsymmL.append( brillouin.BrillouinSegment("Lambda", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "n1",atom_position.Position(0.0,eta,0.0) )) # Lambda line, G to n1
                    self.vsymmL.append( brillouin.BrillouinSegment("C", \
                        "n1",atom_position.Position(1.0,eta,0.0), \
                        "Y",atom_position.Position(1.0,0.0,0.0) )) # C line, n1 to Y
                    self.vsymmL.append( brillouin.BrillouinSegment("l1", \
                        "Y",atom_position.Position(1.0,0.0,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # Y to G
                    self.vsymmL.append( brillouin.BrillouinSegment("l2", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "Z",atom_position.Position(0.0,0.0,0.5) )) # G to Z
                    self.vsymmL.append( brillouin.BrillouinSegment("U", \
                        "Z",atom_position.Position(0.0,0.0,0.5), \
                        "n2",atom_position.Position(0.0,eta,0.5) )) # U line, Z to n2
                    self.vsymmL.append( brillouin.BrillouinSegment("E", \
                        "n2",atom_position.Position(1.0,eta,-0.5), \
                        "M",atom_position.Position(1.0,0.0,-0.5) )) # E line, n2 to M
                    self.vsymmL.append( brillouin.BrillouinSegment("l3", \
                        "M",atom_position.Position(1.0,0.0,-0.5), \
                        "Y",atom_position.Position(1.0,0.0,0.0) )) # M to Y

                else:
                    print("Error: unknown cell type. %s %s %s" % \
                      (self.bravais.shape, self.bravais.center, self.bravais.subtype) )
            else:
                print("Error: unknown cell type. %s %s %s" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        elif self.bravais.shape == "triclinic":
            if self.bravais.center == "simple":
                if self.bravais.subtype == "type1a" or \
                   self.bravais.subtype == "type2a":
                    # Triclinic simple type1a,2a
                    self.vsymmL.append( brillouin.BrillouinSegment("l1", \
                        "X",atom_position.Position(0.5,0.0,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # l line, X to G
                    self.vsymmL.append( brillouin.BrillouinSegment("l2", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "Y",atom_position.Position(0.0,0.5,0.0) )) # l line, G to Y

                    self.vsymmL.append( brillouin.BrillouinSegment("l3", \
                        "L",atom_position.Position(0.5,0.5,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # l line, L to G
                    self.vsymmL.append( brillouin.BrillouinSegment("l4", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "Z",atom_position.Position(0.0,0.0,0.5) )) # l line, G to Z

                    self.vsymmL.append( brillouin.BrillouinSegment("l5", \
                        "N",atom_position.Position(0.5,0.0,0.5), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # l line, N to G
                    self.vsymmL.append( brillouin.BrillouinSegment("l6", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "M",atom_position.Position(0.0,0.5,0.5) )) # l line, G to M

                    self.vsymmL.append( brillouin.BrillouinSegment("l7", \
                        "R",atom_position.Position(0.5,0.5,0.5), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # l line, R to G

                elif self.bravais.subtype == "type1b" or \
                     self.bravais.subtype == "type2b":
                    # Triclinic simple type1b,2b
                    self.vsymmL.append( brillouin.BrillouinSegment("l1", \
                        "X",atom_position.Position(0.0,-0.5,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # l line, X to G
                    self.vsymmL.append( brillouin.BrillouinSegment("l2", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "Y",atom_position.Position(0.5,0.0,0.0) )) # l line, G to Y

                    self.vsymmL.append( brillouin.BrillouinSegment("l3", \
                        "L",atom_position.Position(0.5,-0.5,0.0), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # l line, L to G
                    self.vsymmL.append( brillouin.BrillouinSegment("l4", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "Z",atom_position.Position(-0.5,0.0,0.5) )) # l line, G to Z

                    self.vsymmL.append( brillouin.BrillouinSegment("l5", \
                        #"N",atom_position.Position(-0.5,-0.5,0.5), \
                        "N",atom_position.Position(-0.5,0.5,0.5), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # l line, N to G
                    self.vsymmL.append( brillouin.BrillouinSegment("l6", \
                        "Gamma",atom_position.Position(0.0,0.0,0.0), \
                        "M",atom_position.Position(0.0,0.0,0.5) )) # l line, G to M

                    self.vsymmL.append( brillouin.BrillouinSegment("l7", \
                        "R",atom_position.Position(0.0,-0.5,0.5), \
                        "Gamma",atom_position.Position(0.0,0.0,0.0) )) # l line, R to G

            else:
                print("Error: unknown cell type. %s %s %s" % \
                  (self.bravais.shape, self.bravais.center, self.bravais.subtype) )

        # cut symmetric lines
        s=0
        while s<len(self.vsymmL):
            ks = self.recp_conv[0]*self.vsymmL[s].positions[0] \
               + self.recp_conv[1]*self.vsymmL[s].positions[1] \
               + self.recp_conv[2]*self.vsymmL[s].positions[2]
            ke = self.recp_conv[0]*self.vsymmL[s].positione[0] \
               + self.recp_conv[1]*self.vsymmL[s].positione[1] \
               + self.recp_conv[2]*self.vsymmL[s].positione[2]

            segment = brillouin.BrillouinSegment("","",ks,"",ke)
            cutted = False
            for p in range(len(self.vplane)):
                if not segment.cut(self.vplane[p]):
                    cutted = True
                    break

            if cutted:
                self.vsymmL.pop(s)
                continue

            self.vsymmL[s].positions = atom_position.Position( \
               atom_position.Position_inner_product(self.unit_conv[0],segment.positions), \
               atom_position.Position_inner_product(self.unit_conv[1],segment.positions), \
               atom_position.Position_inner_product(self.unit_conv[2],segment.positions) )

            self.vsymmL[s].positione = atom_position.Position( \
               atom_position.Position_inner_product(self.unit_conv[0],segment.positione), \
               atom_position.Position_inner_product(self.unit_conv[1],segment.positione), \
               atom_position.Position_inner_product(self.unit_conv[2],segment.positione) )
               
            s=s+1

        return True

    def write_kpoints_vasp_line(self, filename, sampling):
        f = open(filename, 'w')
        f.write("band kpath (%s, %s, %s)\n" % \
          (self.bravais.shape, self.bravais.center, self.bravais.subtype))
        f.write("%d\n" % sampling)
        f.write("line\n")
        f.write("reciprocal\n")

        for s in range(len(self.vsymmL)):
            ks = self.recp_conv[0]*self.vsymmL[s].positions[0] \
               + self.recp_conv[1]*self.vsymmL[s].positions[1] \
               + self.recp_conv[2]*self.vsymmL[s].positions[2]

            ke = self.recp_conv[0]*self.vsymmL[s].positione[0] \
               + self.recp_conv[1]*self.vsymmL[s].positione[1] \
               + self.recp_conv[2]*self.vsymmL[s].positione[2]

            ks_prim = atom_position.Position( \
               atom_position.Position_inner_product(self.unit_prim[0],ks), \
               atom_position.Position_inner_product(self.unit_prim[1],ks), \
               atom_position.Position_inner_product(self.unit_prim[2],ks) )
            ke_prim = atom_position.Position( \
               atom_position.Position_inner_product(self.unit_prim[0],ke), \
               atom_position.Position_inner_product(self.unit_prim[1],ke), \
               atom_position.Position_inner_product(self.unit_prim[2],ke) )

            f.write(" %6.3f %6.3f %6.3f   %6.3f %6.3f %6.3f ! %s-%s-%s\n" % \
              (ks_prim[0], ks_prim[1], ks_prim[2], \
               ke_prim[0], ke_prim[1], ke_prim[2], \
               self.vsymmL[s].labelPs, \
               self.vsymmL[s].labelL, \
               self.vsymmL[s].labelPe) )
    
        f.close()

    def write_kpoints(self, filename, total_sampling):
        total_distance = 0.0
        vdistance = []
        vsampling = []
        for s in range(len(self.vsymmL)):
            ks = self.recp_conv[0]*self.vsymmL[s].positions[0] \
               + self.recp_conv[1]*self.vsymmL[s].positions[1] \
               + self.recp_conv[2]*self.vsymmL[s].positions[2]

            ke = self.recp_conv[0]*self.vsymmL[s].positione[0] \
               + self.recp_conv[1]*self.vsymmL[s].positione[1] \
               + self.recp_conv[2]*self.vsymmL[s].positione[2]

            distance = atom_position.Position_distance(ks,ke)
            vdistance.append(distance)
            total_distance += distance

        max_kstep = total_distance/total_sampling

        total_sampling = 0
        for s in range(len(self.vsymmL)):
            sampling = int(math.ceil(vdistance[s]/max_kstep))
            vsampling.append(sampling)
            total_sampling += sampling

            if s+1==len(self.vsymmL) or \
              (s+1< len(self.vsymmL) and \
               not atom_position.Position_match( self.vsymmL[s].positione, \
                                   self.vsymmL[s+1].positions )):
                total_sampling+=1

        f = open( filename, 'w' )
        f.write("band kpath (%s, %s, %s)\n" % \
                    (self.bravais.shape, self.bravais.center, self.bravais.subtype) )
        f.write(" %d\n" % total_sampling )
        f.write("reciprocal\n")

        for s in range(len(self.vsymmL)):
            ks = self.recp_conv[0]*self.vsymmL[s].positions[0] \
               + self.recp_conv[1]*self.vsymmL[s].positions[1] \
               + self.recp_conv[2]*self.vsymmL[s].positions[2]

            ke = self.recp_conv[0]*self.vsymmL[s].positione[0] \
               + self.recp_conv[1]*self.vsymmL[s].positione[1] \
               + self.recp_conv[2]*self.vsymmL[s].positione[2]

            for n in range(vsampling[s]):
                t = float(n)/vsampling[s]
                k = (ke-ks)*t + ks

                k_prim = atom_position.Position( \
                   atom_position.Position_inner_product(self.unit_prim[0],k), \
                   atom_position.Position_inner_product(self.unit_prim[1],k), \
                   atom_position.Position_inner_product(self.unit_prim[2],k) )

                if n==0:
                    f.write("%9.6f %9.6f %9.6f 1.0 ! %d %s %s\n" % \
                                 (k_prim[0],k_prim[1],k_prim[2], \
                                      vsampling[s]-1,self.vsymmL[s].labelPs,self.vsymmL[s].labelL))
                else:
                    f.write("%9.6f %9.6f %9.6f 1.0\n" % (k_prim[0],k_prim[1],k_prim[2]))

            if s+1==len(self.vsymmL) or (s+1<len(self.vsymmL) \
                   and not atom_position.Position_match(self.vsymmL[s].positione,self.vsymmL[s+1].positions)):
                k_prim = atom_position.Position( \
                   atom_position.Position_inner_product(self.unit_prim[0],ke), \
                   atom_position.Position_inner_product(self.unit_prim[1],ke), \
                   atom_position.Position_inner_product(self.unit_prim[2],ke) )

                f.write("%9.6f %9.6f %9.6f 1.0 ! %d %s\n" % \
                            (k_prim[0],k_prim[1],k_prim[2],0,self.vsymmL[s].labelPe))
  
        f.close()

    def write_kpoints_vasp(self, filename, total_sampling):
        kpath = self.get_kpoints(total_sampling)

        f = open(filename,'w')

        f.write("band kpath\n")
        f.write("%d\n" % len(kpath))
        f.write("reciprocal\n")
        for k in kpath:
            if len(k) > 4:
                f.write("%9.6f %9.6f %9.6f %5.2f ! %s\n" % (k[0],k[1],k[2],k[3],k[4]))
            else:
                f.write("%9.6f %9.6f %9.6f %5.2f\n" % (k[0],k[1],k[2],k[3]))
        f.write("\n")  
        f.close()

    def write_kpoints_espresso(self, filename, total_sampling):
        kpath = self.get_kpoints(total_sampling)

        #f = open(filename,'w')
        f = open(filename,'a')

        f.write("\n")
        f.write("K_POINTS crystal\n")
        f.write("%d\n" % len(kpath))
        for k in kpath:
            if len(k) > 4:
                f.write("%9.6f %9.6f %9.6f %5.2f ! %s\n" % (k[0],k[1],k[2],k[3],k[4]))
            else:
                f.write("%9.6f %9.6f %9.6f %5.2f\n" % (k[0],k[1],k[2],k[3]))
        f.write("\n")  
        f.close()

    def write_kpoints_abinit(self, filename, total_sampling):
        kpath = self.get_kpoints(total_sampling)

        #f = open(filename,'w')
        f = open(filename,'a')

        f.write("\n")  
        f.write("nkpt %d\n" % len(kpath))
        f.write("kpt\n")
        for k in kpath:
            if len(k) > 4:
                f.write("%9.6f %9.6f %9.6f # %s\n" % (k[0],k[1],k[2],k[4]))
            else:
                f.write("%9.6f %9.6f %9.6f\n" % (k[0],k[1],k[2]))
        f.write("\n")
        f.close()

    def get_kpoints( self, total_sampling ):
        total_distance = 0.0
        vdistance = []
        vsampling = []
        for s in range(len(self.vsymmL)):
            ks = self.recp_conv[0]*self.vsymmL[s].positions[0] \
               + self.recp_conv[1]*self.vsymmL[s].positions[1] \
               + self.recp_conv[2]*self.vsymmL[s].positions[2]

            ke = self.recp_conv[0]*self.vsymmL[s].positione[0] \
               + self.recp_conv[1]*self.vsymmL[s].positione[1] \
               + self.recp_conv[2]*self.vsymmL[s].positione[2]

            distance = atom_position.Position_distance(ks,ke)
            vdistance.append(distance)
            total_distance += distance

        max_kstep = total_distance/total_sampling

        total_sampling = 0
        for s in range(len(self.vsymmL)):
            sampling = int(math.ceil(vdistance[s]/max_kstep))
            vsampling.append(sampling)
            total_sampling += sampling

            if s+1==len(self.vsymmL) or \
              (s+1< len(self.vsymmL) and \
               not atom_position.Position_match( self.vsymmL[s].positione, \
                                   self.vsymmL[s+1].positions )):
                total_sampling+=1

        kpath = []

        for s in range(len(self.vsymmL)):
            ks = self.recp_conv[0]*self.vsymmL[s].positions[0] \
               + self.recp_conv[1]*self.vsymmL[s].positions[1] \
               + self.recp_conv[2]*self.vsymmL[s].positions[2]

            ke = self.recp_conv[0]*self.vsymmL[s].positione[0] \
               + self.recp_conv[1]*self.vsymmL[s].positione[1] \
               + self.recp_conv[2]*self.vsymmL[s].positione[2]

            for n in range(vsampling[s]):
                t = float(n)/vsampling[s]
                k = (ke-ks)*t + ks

                k_prim = atom_position.Position( \
                   atom_position.Position_inner_product(self.unit_prim[0],k), \
                   atom_position.Position_inner_product(self.unit_prim[1],k), \
                   atom_position.Position_inner_product(self.unit_prim[2],k) )

                if n==0:
                    kpath.append([k_prim[0],k_prim[1],k_prim[2],1.0,self.vsymmL[s].labelPs])
                else:
                    kpath.append([k_prim[0],k_prim[1],k_prim[2],1.0])

            if s+1==len(self.vsymmL) or \
              (s+1<len(self.vsymmL) and \
              not atom_position.Position_match( self.vsymmL[s].positione, \
                                  self.vsymmL[s+1].positions )):
                k_prim = atom_position.Position( \
                   atom_position.Position_inner_product(self.unit_prim[0],ke), \
                   atom_position.Position_inner_product(self.unit_prim[1],ke), \
                   atom_position.Position_inner_product(self.unit_prim[2],ke) )

                kpath.append([k_prim[0],k_prim[1],k_prim[2],1.0,self.vsymmL[s].labelPe])
  
        return kpath

    def write_brillouinzone_gnuplot(self,filename):
        f = open(filename,'w')
        # symmetric points and lines
        f.write("#--- symmetric points and lines\n")
        for s in range(len(self.vsymmL)):
            ks = self.recp_conv[0]*self.vsymmL[s].positions[0] \
                + self.recp_conv[1]*self.vsymmL[s].positions[1] \
                + self.recp_conv[2]*self.vsymmL[s].positions[2]
            ke = self.recp_conv[0]*self.vsymmL[s].positione[0] \
                + self.recp_conv[1]*self.vsymmL[s].positione[1] \
                + self.recp_conv[2]*self.vsymmL[s].positione[2]
            kc = (ks+ke)*0.5

            f.write("set arrow from %9.6f,%9.6f,%9.6f to %9.6f,%9.6f,%9.6f nohead lw 2 lt 2\n" \
                        % (ks[0], ks[1], ks[2], ke[0], ke[1], ke[2]))
            f.write("set label \"%s\" at %9.6f,%9.6f,%9.6f\n" \
                        % (self.vsymmL[s].labelPs, ks[0], ks[1], ks[2]))
            f.write("set label \"%s\" at %9.6f,%9.6f,%9.6f\n" \
                        % (self.vsymmL[s].labelL, kc[0], kc[1], kc[2]))

            if (s+1==len(self.vsymmL)) or \
               (s+1< len(self.vsymmL) and \
                not atom_position.Position_match( self.vsymmL[s].positione, self.vsymmL[s+1].positions)):
                f.write("set label \"%s\" at %9.6f,%9.6f,%9.6f\n" \
                            % (self.vsymmL[s].labelPe,ke[0],ke[1],ke[2]))

        f.write("\n")

        # reciprocal axes
        f.write("#--- reciprocal axis\n")

        label = "a"
        ks = atom_position.Position(0.0,0.0,0.0)
        ke = self.recp_prim[0]
        f.write("set arrow from %9.6f,%9.6f,%9.6f to %9.6f,%9.6f,%9.6f lw 1 lt 3\n" \
                    % (ks[0],ks[1],ks[2],ke[0],ke[1],ke[2]))
        f.write("set label \"%s\" at %+9.6f,%+9.6f,%+9.6f\n" \
                     % (label,ke[0],ke[1],ke[2]))

        label = "b"
        ks = atom_position.Position(0.0,0.0,0.0)
        ke = self.recp_prim[1]
        f.write("set arrow from %9.6f,%9.6f,%9.6f to %9.6f,%9.6f,%9.6f lw 1 lt 3\n" \
                    % (ks[0],ks[1],ks[2],ke[0],ke[1],ke[2]))
        f.write("set label \"%s\" at %9.6f,%9.6f,%9.6f\n" \
                    % (label,ke[0],ke[1],ke[2]))

        label = "c"
        ks = atom_position.Position(0.0,0.0,0.0)
        ke = self.recp_prim[2]
        f.write("set arrow from %9.6f,%9.6f,%9.6f to %9.6f,%9.6f,%9.6f lw 1 lt 3\n" \
                    % (ks[0],ks[1],ks[2],ke[0],ke[1],ke[2]))
        f.write("set label \"%s\" at %9.6f,%9.6f,%9.6f\n" \
                    % (label,ke[0],ke[1],ke[2]))

        f.write("\n")

        for i in range(len(self.vfacet)):
            f.write("#--- facet %d\n" % i)
            self.vfacet[i].output(f)
        f.write("\n")

        f.write("set nokey\n")
        f.write("unset tics\n")
        f.write("set ticslevel 0\n")
        f.write("unset border\n")
        f.write("set size square\n")
        f.write("set view 60,110\n")
        f.write("\n")

        f.write("set title \"%s %s %s\"\n" \
                     % (self.bravais.shape,self.bravais.center,self.bravais.subtype))

        scale=0.0
        for i in [0,1,2]:
            v = abs(self.recp_prim[i][0])
            if scale < v: scale = v
            v = abs(self.recp_prim[i][1])
            if scale < v: scale = v
            v = abs(self.recp_prim[i][2])
            if scale < v: scale = v

        scale *= 0.75

        f.write("set xrange [%9.6f:%9.6f]\n" % (-scale,scale))
        f.write("set yrange [%9.6f:%9.6f]\n" % (-scale,scale))
        f.write("set zrange [%9.6f:%9.6f]\n" % (-scale,scale))

        os.system("touch gamma.dat")
        f.write("splot \"gamma.dat\"\n")

        #f.write("set terminal png\n")
        #f.write("set output \"tmp.png\"\n")
        #f.write("replot\n")

        f.write("pause -1 \"hit return key\"\n")

        f.close()

    def write_brillouinzone_jmol(self,filename):
        f = open(filename,'w')
        # symmetric points and lines
        f.write("#--- symmetric points and lines\n")
        for s in range(len(self.vsymmL)):
            ks = self.recp_conv[0]*self.vsymmL[s].positions[0] \
                + self.recp_conv[1]*self.vsymmL[s].positions[1] \
                + self.recp_conv[2]*self.vsymmL[s].positions[2]
            ke = self.recp_conv[0]*self.vsymmL[s].positione[0] \
                + self.recp_conv[1]*self.vsymmL[s].positione[1] \
                + self.recp_conv[2]*self.vsymmL[s].positione[2]
            kc = (ks+ke)*0.5

            f.write("draw l%d \"\" {%9.6f %9.6f %9.6f} {%9.6f %9.6f %9.6f};\n" \
                        % (s,ks[0],ks[1],ks[2],ke[0],ke[1],ke[2]) )
            f.write("draw p0%d \"%s\" {%9.6f %9.6f %9.6f}\n" \
                        % (s,self.vsymmL[s].labelPs,ks[0], ks[1], ks[2]) )

            if (s+1==len(self.vsymmL)) or \
               (s+1< len(self.vsymmL) and \
                not atom_position.Position_match( self.vsymmL[s].positione, self.vsymmL[s+1].positions )):

                f.write("draw p1%d \"%s\" {%9.6f %9.6f %9.6f}\n" \
                            % (s,self.vsymmL[s].labelPe,ke[0], ke[1],ke[2]))

        f.write("\n")

        # reciprocal axes
        f.write("#--- reciprocal axis\n")

        ks = atom_position.Position(0.0,0.0,0.0)
        ke = self.recp_prim[0]

        f.write("draw a \"a\" {%9.6f %9.6f %9.6f} {%9.6f %9.6f %9.6f} color red;\n" \
                    % (ke[0],ke[1],ke[2],ks[0],ks[1],ks[2]))

        ks = atom_position.Position(0.0,0.0,0.0)
        ke = self.recp_prim[1]

        f.write("draw b \"b\" {%9.6f %9.6f %9.6f} {%9.6f %9.6f %9.6f} color green;\n" \
                    % (ke[0],ke[1],ke[2],ks[0],ks[1],ks[2]))

        ks = atom_position.Position(0.0,0.0,0.0)
        ke = self.recp_prim[2]

        f.write("draw c \"c\" {%9.6f %9.6f %9.6f} {%9.6f %9.6f %9.6f} color blue;\n" \
                    % (ke[0],ke[1],ke[2],ks[0],ks[1],ks[2]))

        f.write("\n")

        for i in range(len(self.vfacet)):
            f.write("#--- begin facet %d\n" % i)
            self.vfacet[i].output2(f,i)

        f.write("\n")

        scale=0.0
        for i in [0,1,2]:
            v = abs(self.recp_prim[i][0])
            if scale < v: scale = v
            v = abs(self.recp_prim[i][1])
            if scale < v: scale = v
            v = abs(self.recp_prim[i][2])
            if scale < v: scale = v

        scale *= 0.75
        f.write("zoom %12.2f;\n" % (10.0/scale/2*100) )

        f.close()

def make_symop(rot, trans):
    '''make symop symbols from rotation matrix and translation vectors
       rot: rotational matrix in numpy 3x3 array
       trans: translational vector
       return string 
    '''
    symbols = []
    for i in range(3):
        symop = ""
        if rot[i,0] == 1:
            symop += "x"
        elif rot[i,0] == -1:
            symop += "-x"
        elif rot[i,0] == 0:
            pass
        else:
            print("error in rotation matrix")
        if rot[i,1] == 1:
            if symop:
                symop += "+y"
            else:
                symop += "y"
        elif rot[i,1] == -1:
            symop += "-y"
        elif rot[i,1] == 0:
            pass
        else:
            print("error in rotation matrix")
        if rot[i,2] == 1:
            if symop:
                symop += "+z"
            else:
                symop += "z"
        elif rot[i,2] == -1:
            symop += "-z"
        elif rot[i,2] == 0:
            pass
        else:
            print("error in rotation matrix")
        (shift, dummy) = math.modf(trans[i])
        if shift < 0:
            shift += 1.0
        if abs(shift) > 1.0e-4 and abs(shift-1.0) > 1.0e-4:
            if shift > 0:
                symop += "+"
            # symop += str(trans[i])
            symop += float2rational(trans[i])
        symbols.append(symop)
    return ",".join(symbols)

def pearsonsymbol(spgnum):
    spglist = spgtable.SpaceGroupTable()
    spg = spglist.findByIntTableNumber(spgnum)
    family = None
    if spgnum < 3:
        family = "a"
    elif spgnum < 16:
        family = "m"
    elif spgnum < 75:
        family = "o"
    elif spgnum < 143:
        family = "t"
    elif spgnum < 195:
        family = "h"
    else:
        family = "c"
    centring = spg.name[0]
    if centring in "ABC":
        centring = "S"
    return family+centring

def roundint(x):
    return int(math.floor(x+0.5))

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
        # print(x, d, w, n, w-n)
        if abs(w-n) < tol:
            #s = "{:d}/{:d}".format(n, d)
            s = "%d/%d" % (n, d)
            # print("float2rational:", x, s)
            return s
    sys.stderr.write("can't convert to rational {}\n".format(x))
    sys.exit(1)

def usage( argv ):
    print("usage:")

    print(" case A: CIF to VASP(POSCAR,KPOINTS) in primitive cell")
    print("   $ %s -f input.cif POSCAR KPOINTS plot_bz.plt plot_bz.spt" % argv[0] )
    print("        input.cif is the input CIF file.")
    print("        POSCAR is the output CAR file in primitive cell.")
    print("        KPOINTS is the output KPOINTS file.")
    print("        plot_bz.plt is the output gnuplot script (brillouin zone).")
    print("        plot_bz.spt is the output jmol script (brillouin zone).")
    print("\n")

    print(" case B: VASP(POSCAR) to CIF(P1) and VASP(KPOINTS)")
    print("   $ %s -b2 POSCAR output.cif KPOINTS plot_bz.plt plot_bz.spt" % argv[0] )
    print("        POSCAR is the input CAR file.")
    print("        output.cif is the output CIF file (P1).")
    print("        KPOINTS is the output KPOINTS file.")
    print("        plot_bz.plt is the output gnuplot script. (brillouin zone)")
    print("        plot_bz.spt is the output jmol script. (brillouin zone)")
    print("\n")

    print(" case C: VASP(POSCAR) and CIF to another CIF")
    print("   $ %s -b POSCAR input.cif output.cif" % argv[0] )
    print("        POSCAR is the input CAR file.")
    print("        input.cif is the input CIF file.")
    print("        output.cif is the output CIF file.")
    print("\n")

    print(" case D: CIF to VASP(POSCAR,KPOINTS) in conventional cell")
    print("   $ %s -f2 input.cif POSCAR KPOINTS plot_bz.plt POSCAR2  plot_bz.spt" % argv[0] )
    print("        input.cif is the input CIF file.")
    print("        POSCAR is the output CAR file in primitive cell.")
    print("        KPOINTS is the output KPOINTS file.")
    print("        plot_bz.plt is the output gnuplot script (brillouin zone).")
    print("        POSCAR2 is the output CAR file in conventional cell.")
    print("        plot_bz.spt is the output jmol script (brillouin zone).")
    print("\n")

    print(" case E: CIF to BZ plot file")
    print("   $ %s -bz input.cif plot_bz.plt plot_bz.spt" % argv[0] )
    print("        input.cif is the input CIF file.")
    print("        plot_bz.plt is the output gnuplot script (brillouin zone).")
    print("        plot_bz.spt is the output jmol script (brillouin zone).")

def csa2():
    arg = argparse.ArgumentParser(description="crystal analyzer")
    arg.add_argument("-c","--cif",action="store",dest="cif")
    arg.add_argument("--poscar",action="store",dest="poscar")
    arg.add_argument("-p","--program",action="store",dest="program",default="qe")
    arg.add_argument("-o","--output",action="store",dest="fout")
    arg.add_argument("--ocif",action="store",dest="foutcif")
    arg.add_argument("--bz_gnuplot",action="store",dest="bz_gnuplot")
    arg.add_argument("--bz_jmol",action="store",dest="bz_jmol")
    arg.add_argument("--nkpath",action="store",dest="nkpath",default=100)
    args = arg.parse_args()
    cry = Crystal()
    print
    print("crystal analyzer")
    print("ver.",cry.version)
    print("")
    print("[crystal structure analysis]")

    if args.cif != None:
        print("cif:",args.cif)
        cry.read_cif(args.cif)
        cry.check_cif()

    cry.construct_structure()

    if args.poscar != None:
        print("poscar:",args.poscar)
        cry.read_structure_vasp(args.poscar)        
        cry.reconstruct_structure()

    cry.print_log()

    print("program:",args.program)

    total_sampling = 100
    if args.nkpath != None:
        total_sampling = int(args.nkpath)

    if args.program=="vasp":
        print("output: POSCAR, KPOINTS")
        cry.write_structure_vasp("POSCAR")
        cry.write_kpoints_vasp("KPOINTS",total_sampling)
        cry.write_kpoints_vasp_line("KPOINTS.line",10)
    elif args.program=="qe":
        if args.fout != None:
            print("output:",args.fout)
            cry.write_structure_espresso(args.fout)
            cry.write_kpoints_espresso(args.fout,total_sampling)
    elif args.program=="abinit":
        if args.fout != None:
            print("output:",args.fout)
            cry.write_structure_abinit(args.fout)
            cry.write_kpoints_abinit(args.fout,total_sampling)

    if args.foutcif != None:
        print("output cif:",args.foutcif)
        cry.write_cif(args.foutcif)

    if args.bz_gnuplot != None:
        print("output bz (gnuplot):",args.bz_gnuplot)
        cry.write_brillouinzone_gnuplot(args.bz_gnuplot)
    if args.bz_jmol != None:
        print("output bz (jmol):",args.bz_jmol)
        cry.write_brillouinzone_jmol(args.bz_jmol)

    cry.write_cif2("tmp2.cif")
#------------------------------------------------------------- 
if __name__ == "__main__":
    csa2()
