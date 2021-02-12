import sys
import copy
from. import position 

class BrillouinPlane():
    def __init__( self, normal=position.Position() ):
        self.normal = normal
        self.distance = (-1.0)*position.Position_inner_product(self.normal,self.normal)

    def equation( self, atomic_position):
        return position.Position_inner_product(self.normal,atomic_position) + self.distance

class BrillouinSegment():
    def __init__( self, labelL="", \
            labelPs="", positions=position.Position(), \
            labelPe="", positione=position.Position() ):
        self.labelL    = labelL
        self.labelPs   = labelPs
        self.positions = positions
        self.labelPe   = labelPe
        self.positione = positione
        self.positionc = position.Position()
        self.cutted    = False

    def cut( self, plane ):
        vs = plane.equation(self.positions)
        ve = plane.equation(self.positione)

        eps = 1.0e-8

        # the both point is on the plane
        if abs(vs)<=eps and abs(ve)<=eps:
            self.cutted = False
            return True

        # the end point is on the plane, the segment is inside
        if vs< 0.0 and abs(ve)<=eps:
            self.positionc = self.positione
            self.cutted = True
            return True

        # the end point is on the plane, the segment is outside
        if vs> 0.0 and abs(ve)<=eps:
            self.positionc = self.positione
            self.cutted = True
            return False

        # the start point is on the plane, the segment is outside
        if abs(vs)<=eps and ve>0.0:
            self.positionc = self.positions
            self.cutted = True
            return False

        # the start point is on the plane, the segment is inside
        if abs(vs)<=eps and ve<0.0:
            self.positionc = self.positions
            self.cutted = True
            return True

        # the segment is inside the plane
        if vs< 0.0 and ve<0.0:
            self.cutted = False
            return True

        # the segment crosses the plane. vs outside, ve inside.
        if vs> 0.0 and ve<0.0:
            t  = (vs)/(vs-ve)
            self.positionc = self.positione*t + self.positions*(1.0-t)
            self.positions = self.positionc
            self.cutted = True
            return True

        # the segment crosses the plane.
        if vs< 0.0 and ve>0.0:
            t = (ve)/(ve-vs)
            self.positionc = self.positions*t + self.positione*(1.0-t)
            self.positione = self.positionc
            self.cutted = True
            return True

        # the segment is outside the plane
        if vs> 0.0 and ve>0.0:
            self.cutted = False
            return False

        print("# Abnormal error: BrillouinSegment::cut.")
        sys.exit(0)
        return False

    def output( self, fd ):
        fd.write( "set arrow from %+6.3f,%+6.3f,%+6.3f to %+6.3f,%+6.3f,%+6.3f nohead lw 1 lt 1\n" % \
          ( self.positions[0], self.positions[1], self.positions[2], \
            self.positione[0], self.positione[1], self.positione[2] ) )

    def output2( self, fd, f, s ):
        fd.write( "draw f%d%d {%6.3f %6.3f %6.3f} {%6.3f %6.3f %6.3f} color white;\n" % \
          ( f, s, \
            self.positions[0], self.positions[1], self.positions[2], \
            self.positione[0], self.positione[1], self.positione[2] ) )

def BrillouinSegment_match( segment1, segment2 ):
    return \
    ( position.Position_match(segment1.positions,segment2.positions) and \
      position.Position_match(segment1.positione,segment2.positione) ) or \
    ( position.Position_match(segment1.positions,segment2.positione) and \
      position.Position_match(segment1.positione,segment2.positions) )

def BrillouinSegment_find( vsegment, segment ):
    for n in range(len(vsegment)):
        if BrillouinSegment_match(vsegment[n],segment):
            return True
    return False

class BrillouinFacet():
    def __init__( self, \
      position0=position.Position(), \
      position1=position.Position(), \
      position2=position.Position(), \
      position3=position.Position() ):
        self.vsegment = []
        self.vsegment.append( BrillouinSegment("","",position0,"",position1) )
        self.vsegment.append( BrillouinSegment("","",position1,"",position2) )
        self.vsegment.append( BrillouinSegment("","",position2,"",position3) )
        self.vsegment.append( BrillouinSegment("","",position3,"",position0) )
        self.segmentc = BrillouinSegment("","",position.Position(),"",position.Position())
        self.cutted = False

    def cut( self, plane ):
        cutted_size = 0

        s=0
        while s<len(self.vsegment):
            inside = self.vsegment[s].cut(plane)

            if self.vsegment[s].cutted:
                if cutted_size == 0:
                    self.segmentc.positions = self.vsegment[s].positionc
                    cutted_size+=1

                elif cutted_size == 1:
                    if not position.Position_match(self.segmentc.positions,self.vsegment[s].positionc):
                        self.segmentc.positione = self.vsegment[s].positionc
                        cutted_size+=1

            if not inside:
                self.vsegment.pop(s)
                continue
            s+=1

        if cutted_size == 2:
            self.vsegment.append(copy.deepcopy(self.segmentc))
            self.cutted = True
        else:
            self.cutted = False

        if len(self.vsegment)<3:
            return False

        if self.cutted:
            self.order()

        return True

    def order( self ):
        s=0
        while s<len(self.vsegment):
            t=s+1
            while t<len(self.vsegment):
                if BrillouinSegment_match( self.vsegment[s], self.vsegment[t] ):
                    self.vsegment.pop( t )
                    continue
                t+=1
            s+=1

        for s in range(len(self.vsegment)):
            if position.Position_match( self.vsegment[s].positions, self.vsegment[s].positione ):
                print("# Abnormal error: BrillouinFacet::order. s==e.")
                sys.exit(0)

        for s in range(len(self.vsegment)):
            for t in range(s+1,len(self.vsegment)):
                if position.Position_match( self.vsegment[s].positione, self.vsegment[t].positions ):
                    if s+1 != t:
                        self.vsegment[s+1], self.vsegment[t] = \
                        self.vsegment[t], self.vsegment[s+1]
                    break
                if position.Position_match( self.vsegment[s].positione, self.vsegment[t].positione ):
                    self.vsegment[t].positions, self.vsegment[t].positione = \
                    self.vsegment[t].positione, self.vsegment[t].positions

                    if s+1 != t:
                        self.vsegment[s+1], self.vsegment[t] = \
                        self.vsegment[t], self.vsegment[s+1]
                    break

        s=0
        e=len(self.vsegment)-1
        if not position.Position_match( self.vsegment[s].positions, self.vsegment[e].positione ):
            print("# Abnormal error: BrillouinFacet::order. not loop.")
            for s in range(len(self.vsegment)):
                self.vsegment[s].output(sys.stdout)
            sys.exit(0)

    def output( self, fd ):
        for s in range(len(self.vsegment)):
            self.vsegment[s].output(fd)

    def output2( self, fd, f ):
        for s in range(len(self.vsegment)):
            self.vsegment[s].output2(fd,f,s)
