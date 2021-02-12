from . import position as pos 

class Atom():
    def __init__(self):
        self.label = ''
        self.position = [0.0,0.0,0.0]
        self.occupancy = 0.0
        self.multiplicity = 0
        self.wyckoff = ''
        self.uiso = '0.0'
        self.element = ''
        self.site_index = 0
        self.operation_index = 0

def Atom_match(atom1, atom2):
    eps = 2.0e-4
    return pos.Position_match2( atom1.position, atom2.position,eps)

def Atom_find(vatom, atom):
    if vatom==[]:
        return False
    else:
        for a in vatom:
            if Atom_match(a,atom):
                return True
            else:
                pass
        return False
