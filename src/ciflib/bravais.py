import math

class Bravais():
    def __init__(self):
        self.shape = ''
        self.center = ''
        self.subtype = ''

    def check(self, spacegroup, cell):
        self.center = spacegroup.center
        eps = 1.0e-5

        if spacegroup.shape == 'cubic':
            if(cell.length_a == cell.length_b and \
               cell.length_b == cell.length_c and \
               cell.length_c == cell.length_a and \
               abs(cell.angle_alpha-90.0)<eps and \
               abs(cell.angle_beta -90.0)<eps and \
               abs(cell.angle_gamma-90.0)<eps):
                self.shape = 'cubic'
                self.subtype = 'none'
                return True
            else:
                return False

        elif spacegroup.shape == 'tetragonal':
            if(cell.length_a == cell.length_b and \
               abs(cell.angle_alpha-90.0)<eps and \
               abs(cell.angle_beta -90.0)<eps and \
               abs(cell.angle_gamma-90.0)<eps):
                self.shape = 'tetragonal'
                if(self.center == 'simple'):
                    self.subtype = 'none'
                elif self.center == 'body':
                    if(cell.length_a < cell.length_c):
                        self.subtype = 'type1'
                    else:
                        self.subtype = 'type2'
                return True
            else:
                return False

        elif spacegroup.shape == 'trigonal':
            if(cell.length_a == cell.length_b and \
               abs(cell.angle_alpha-90.0)<eps and \
               abs(cell.angle_beta -90.0)<eps and \
               abs(cell.angle_gamma-120.0)<eps):
                # as hexagonal conventional cell
                self.shape = 'trigonal'
                if(cell.length_a*math.sqrt(3)<cell.length_c*math.sqrt(2)):
                    self.subtype = 'type1'
                else:
                    self.subtype = 'type2'
                return True

            elif(cell.length_a == cell.length_b and \
                 cell.length_b == cell.length_c and \
                 cell.length_c == cell.length_a and \
                 cell.angle_gamma == cell.angle_alpha and \
                 cell.angle_alpha == cell.angle_beta):
                # as trigonal conventional cell
                self.shape = 'trigonal'
                if(cell.angle_alpha<90.0):
                    self.subtype = 'type1'
                else:
                    self.subtype = 'type2'
                return True
            else:
                return False

        elif(spacegroup.shape == 'hexagonal'):
            if(cell.length_a == cell.length_b and \
               abs(cell.angle_alpha-90.0)<eps and \
               abs(cell.angle_beta -90.0)<eps and \
               abs(cell.angle_gamma-120.0)<eps):
                self.shape = 'hexagonal'
                self.subtype = 'none'
                return True
            else:
                return False

        elif(spacegroup.shape == 'orthorhombic'):
            if(abs(cell.angle_alpha-90.0)<eps and \
               abs(cell.angle_beta -90.0)<eps and \
               abs(cell.angle_gamma-90.0)<eps):
                self.shape = 'orthorhombic'

                if(self.center == 'simple'):
                    self.subtype = 'none'
                elif(self.center == 'face'):
                    ka2 = 1.0/(cell.length_a*cell.length_a)
                    kb2 = 1.0/(cell.length_b*cell.length_b)
                    kc2 = 1.0/(cell.length_c*cell.length_c)
                    if ka2 < kb2 + kc2 and \
                       kb2 < ka2 + kc2 and \
                       kc2 < ka2 + kb2:
                        self.subtype = 'type1'
                    elif ka2 > kb2 + kc2:
                        self.subtype = 'type2'
                    elif kb2 > ka2 + kc2:
                        self.subtype = 'type3'
                    elif kc2 > ka2 + kb2:
                        self.subtype = 'type4'
                elif(self.center == 'body'):
                    if(cell.length_a>cell.length_c and \
                         cell.length_a>cell.length_b):
                        self.subtype = 'type1'
                    elif(cell.length_b>cell.length_a and \
                         cell.length_b>cell.length_c):
                        self.subtype = 'type2'
                    elif(cell.length_c>cell.length_a and \
                         cell.length_c>cell.length_b):
                        self.subtype = 'type3'
                elif(self.center == 'base'):
                    # suppose C base
                    if(cell.length_a<cell.length_b):
                        self.subtype = 'type1'
                    else:
                        self.subtype = 'type2'
                return True
            else:
                return False

        elif(spacegroup.shape == 'monoclinic'):
            if(abs(cell.angle_gamma-90.0)<eps and \
               abs(cell.angle_alpha-90.0)<eps):
                self.shape = 'monoclinic'
                if(self.center == 'simple'):
                    self.subtype = 'none'
                elif(self.center == 'base'): # AB plane is the base
                    if(cell.length_a<cell.length_b):
                        self.subtype = 'type1'
                    else:
                        self.subtype = 'type2'
                return True
            else:
                return False

        elif(spacegroup.shape == 'triclinic'):
            self.shape = 'triclinic'
            self.subtype = 'none' # can not determine here
            return True

        self.shape = 'none'
        self.subtype = 'none'

        return False
