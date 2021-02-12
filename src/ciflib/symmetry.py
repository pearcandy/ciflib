#!/usr/bin/env python
import re
from . import position as pos

class SymmetryOperation():
    def __init__(self,codes):
        self.codes = codes
        self.transM = [[0,0,0],[0,0,0],[0,0,0]]
        self.transV = [0,0,0]

        m = re.match('([-+xyz]+)([/0-9]*),[ ]*([-+xyz]+)([/0-9]*),[ ]*([-+xyz]+)([/0-9]*)',codes)
        if(m):
            axis = [m.group(1),m.group(3),m.group(5)]
            frac = [m.group(2),m.group(4),m.group(6)]
        else:
            m = re.match('([/0-9]*)([-+xyz]+),[ ]*([/0-9]*)([-+xyz]+),[ ]*([/0-9]*)([-+xyz]+)',codes)
            axis = [m.group(2),m.group(4),m.group(6)]
            frac = [m.group(1),m.group(3),m.group(5)]

        for i in range(3):
            if(axis[i].find('-x')>=0):
                self.transM[0][i] = -1.0
            elif(axis[i].find('x')>=0):
                self.transM[0][i] = +1.0
            else:
                self.transM[0][i] =  0.0

            if(axis[i].find('-y')>=0):
                self.transM[1][i] = -1.0
            elif(axis[i].find('y')>=0):
                self.transM[1][i] = +1.0
            else:
                self.transM[1][i] =  0.0

            if(axis[i].find('-z')>=0):
                self.transM[2][i] = -1.0
            elif(axis[i].find('z')>=0):
                self.transM[2][i] = +1.0
            else:
                self.transM[2][i] =  0.0

            if(frac[i] != ''):
                m = re.match('(\d+)/(\d+)',frac[i])
                if(m):
                    self.transV[i] = float(m.group(1))/float(m.group(2))
                else:
                    self.transV[i] = 0.0

    def operate(self, p_org):
        p_new = pos.Position_translate(self.transM, p_org)

        p_new[0] += self.transV[0]
        p_new[1] += self.transV[1]
        p_new[2] += self.transV[2]

        pos.Position_fold(p_new)

        return p_new
