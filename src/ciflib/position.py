import math
import numpy

def Position(x=0.0,y=0.0,z=0.0):
    return numpy.array([x,y,z])

def Position_match(p1, p2, eps=1.0e-8):
    return \
      abs(p1[0]-p2[0])<eps and \
      abs(p1[1]-p2[1])<eps and \
      abs(p1[2]-p2[2])<eps

def Position_match2(p1, p2, eps=1.0e-8):
    dx = abs(p1[0]-p2[0])
    if (dx>eps and abs(dx-1.0)>eps):
        return False
    dy = abs(p1[1]-p2[1])
    if (dy>eps and abs(dy-1.0)>eps):
        return False
    dz = abs(p1[2]-p2[2])
    if (dz>eps and abs(dz-1.0)>eps):
        return False
    return True

def Position_outer_product(p1, p2):
    x = p1[1]*p2[2] - p2[1]*p1[2]
    y = p1[2]*p2[0] - p2[2]*p1[0]
    z = p1[0]*p2[1] - p2[0]*p1[1]
    return Position(x,y,z)

def Position_inner_product(p1, p2):
    return p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2]

def Position_distance(p1, p2):
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    dz = p1[2] - p2[2]
    return math.sqrt(dx*dx+dy*dy+dz*dz)

def Position_length(p):
    dx = p[0]
    dy = p[1]
    dz = p[2]
    return math.sqrt(dx*dx+dy*dy+dz*dz)

def Position_angle(p1, p2):
    len1 = Position_length(p1)
    len2 = Position_length(p2)
    prod = Position_inner_product(p1,p2)
    return math.degrees(math.acos(prod/(len1*len2)))

def Position_reciprocal(unit):
    recp = [Position(), Position(), Position()]
    recp[0] = Position_outer_product(unit[1], unit[2])
    recp[1] = Position_outer_product(unit[2], unit[0])
    recp[2] = Position_outer_product(unit[0], unit[1])

    Vunit = Position_inner_product(recp[0], unit[0])

    recp[0] = recp[0]*(1.0/Vunit) # 2.0*M_PI not multiplied
    recp[1] = recp[1]*(1.0/Vunit) # 2.0*M_PI not multiplied
    recp[2] = recp[2]*(1.0/Vunit) # 2.0*M_PI not multiplied

    return recp

def Position_translate(vecs, p_org):
    x = vecs[0][0]*p_org[0] + vecs[1][0]*p_org[1] + vecs[2][0]*p_org[2]
    y = vecs[0][1]*p_org[0] + vecs[1][1]*p_org[1] + vecs[2][1]*p_org[2]
    z = vecs[0][2]*p_org[0] + vecs[1][2]*p_org[1] + vecs[2][2]*p_org[2]
    return Position(x,y,z)

def Position_fold(p):
    while p[0]< 0.0: p[0]+=1.0
    while p[0]>=1.0: p[0]-=1.0
    while p[1]< 0.0: p[1]+=1.0
    while p[1]>=1.0: p[1]-=1.0
    while p[2]< 0.0: p[2]+=1.0
    while p[2]>=1.0: p[2]-=1.0
