"""
Copyright (C) 2018 National Institute for Materials Science (NIMS)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE."""


class SpaceGroup():
    def __init__(self,
                 hallnumber=0,
                 number=0,
                 name='',
                 nameHM='',
                 nameFull='',
                 nameHall='',
                 axis_code='',
                 pointgroup='',
                 shape='',
                 center=''):
        self.hallnumber = hallnumber
        self.number = number
        self.name = name
        self.nameHM = nameHM
        self.nameFull = nameFull
        self.nameHall = nameHall
        self.axis_code = axis_code
        self.pointgroup = pointgroup
        self.shape = shape
        self.center = center


class SpaceGroupTable():
    def __init__(self):
        self.spacegroup = []
        self.spacegroup.append(
            SpaceGroup(1, 1, 'P1', 'P 1', 'P 1', 'P 1', '', '1', 'triclinic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(2, 2, 'P-1', 'P -1', 'P -1', '-P 1', '', '-1',
                       'triclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(3, 3, 'P2', 'P 2 = P 1 2 1', 'P 1 2 1', 'P 2y', 'b',
                       '2', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(4, 3, 'P2', 'P 2 = P 1 1 2', 'P 1 1 2', 'P 2', 'c', '2',
                       'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(5, 3, 'P2', 'P 2 = P 2 1 1', 'P 2 1 1', 'P 2x', 'a',
                       '2', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(6, 4, 'P2_1', 'P 2_1 = P 1 2_1 1', 'P 1 2_1 1', 'P 2yb',
                       'b', '2', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(7, 4, 'P2_1', 'P 2_1 = P 1 1 2_1', 'P 1 1 2_1', 'P 2c',
                       'c', '2', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(8, 4, 'P2_1', 'P 2_1 = P 2_1 1 1', 'P 2_1 1 1', 'P 2xa',
                       'a', '2', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(9, 5, 'C2', 'C 2 = C 1 2 1', 'C 1 2 1', 'C 2y', 'b1',
                       '2', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(10, 5, 'C2', 'C 2 = A 1 2 1', 'A 1 2 1', 'A 2y', 'b2',
                       '2', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(11, 5, 'C2', 'C 2 = I 1 2 1', 'I 1 2 1', 'I 2y', 'b3',
                       '2', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(12, 5, 'C2', 'C 2 = A 1 1 2', 'A 1 1 2', 'A 2', 'c1',
                       '2', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(13, 5, 'C2', 'C 2 = B 1 1 2 = B 2', 'B 1 1 2', 'B 2',
                       'c2', '2', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(14, 5, 'C2', 'C 2 = I 1 1 2', 'I 1 1 2', 'I 2', 'c3',
                       '2', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(15, 5, 'C2', 'C 2 = B 2 1 1', 'B 2 1 1', 'B 2x', 'a1',
                       '2', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(16, 5, 'C2', 'C 2 = C 2 1 1', 'C 2 1 1', 'C 2x', 'a2',
                       '2', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(17, 5, 'C2', 'C 2 = I 2 1 1', 'I 2 1 1', 'I 2x', 'a3',
                       '2', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(18, 6, 'Pm', 'P m = P 1 m 1', 'P 1 m 1', 'P -2y', 'b',
                       'm', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(19, 6, 'Pm', 'P m = P 1 1 m', 'P 1 1 m', 'P -2', 'c',
                       'm', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(20, 6, 'Pm', 'P m = P m 1 1', 'P m 1 1', 'P -2x', 'a',
                       'm', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(21, 7, 'Pc', 'P c = P 1 c 1', 'P 1 c 1', 'P -2yc', 'b1',
                       'm', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(22, 7, 'Pc', 'P c = P 1 n 1', 'P 1 n 1', 'P -2yac',
                       'b2', 'm', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(23, 7, 'Pc', 'P c = P 1 a 1', 'P 1 a 1', 'P -2ya', 'b3',
                       'm', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(24, 7, 'Pc', 'P c = P 1 1 a', 'P 1 1 a', 'P -2a', 'c1',
                       'm', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(25, 7, 'Pc', 'P c = P 1 1 n', 'P 1 1 n', 'P -2ab', 'c2',
                       'm', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(26, 7, 'Pc', 'P c = P 1 1 b = P b', 'P 1 1 b', 'P -2b',
                       'c3', 'm', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(27, 7, 'Pc', 'P c = P b 1 1', 'P b 1 1', 'P -2xb', 'a1',
                       'm', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(28, 7, 'Pc', 'P c = P n 1 1', 'P n 1 1', 'P -2xbc',
                       'a2', 'm', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(29, 7, 'Pc', 'P c = P c 1 1', 'P c 1 1', 'P -2xc', 'a3',
                       'm', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(30, 8, 'Cm', 'C m = C 1 m 1', 'C 1 m 1', 'C -2y', 'b1',
                       'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(31, 8, 'Cm', 'C m = A 1 m 1', 'A 1 m 1', 'A -2y', 'b2',
                       'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(32, 8, 'Cm', 'C m = I 1 m 1', 'I 1 m 1', 'I -2y', 'b3',
                       'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(33, 8, 'Cm', 'C m = A 1 1 m', 'A 1 1 m', 'A -2', 'c1',
                       'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(34, 8, 'Cm', 'C m = B 1 1 m = B m', 'B 1 1 m', 'B -2',
                       'c2', 'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(35, 8, 'Cm', 'C m = I 1 1 m', 'I 1 1 m', 'I -2', 'c3',
                       'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(36, 8, 'Cm', 'C m = B m 1 1', 'B m 1 1', 'B -2x', 'a1',
                       'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(37, 8, 'Cm', 'C m = C m 1 1', 'C m 1 1', 'C -2x', 'a2',
                       'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(38, 8, 'Cm', 'C m = I m 1 1', 'I m 1 1', 'I -2x', 'a3',
                       'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(39, 9, 'Cc', 'C c = C 1 c 1', 'C 1 c 1', 'C -2yc', 'b1',
                       'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(40, 9, 'Cc', 'C c = A 1 n 1', 'A 1 n 1', 'A -2yac',
                       'b2', 'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(41, 9, 'Cc', 'C c = I 1 a 1', 'I 1 a 1', 'I -2ya', 'b3',
                       'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(42, 9, 'Cc', 'C c = A 1 a 1', 'A 1 a 1', 'A -2ya',
                       '-b1', 'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(43, 9, 'Cc', 'C c = C 1 n 1', 'C 1 n 1', 'C -2ybc',
                       '-b2', 'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(44, 9, 'Cc', 'C c = I 1 c 1', 'I 1 c 1', 'I -2yc',
                       '-b3', 'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(45, 9, 'Cc', 'C c = A 1 1 a', 'A 1 1 a', 'A -2a', 'c1',
                       'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(46, 9, 'Cc', 'C c = B 1 1 n', 'B 1 1 n', 'B -2bc', 'c2',
                       'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(47, 9, 'Cc', 'C c = I 1 1 b', 'I 1 1 b', 'I -2b', 'c3',
                       'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(48, 9, 'Cc', 'C c = B 1 1 b = B b', 'B 1 1 b', 'B -2b',
                       '-c1', 'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(49, 9, 'Cc', 'C c = A 1 1 n', 'A 1 1 n', 'A -2ac',
                       '-c2', 'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(50, 9, 'Cc', 'C c = I 1 1 a', 'I 1 1 a', 'I -2a', '-c3',
                       'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(51, 9, 'Cc', 'C c = B b 1 1', 'B b 1 1', 'B -2xb', 'a1',
                       'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(52, 9, 'Cc', 'C c = C n 1 1', 'C n 1 1', 'C -2xbc',
                       'a2', 'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(53, 9, 'Cc', 'C c = I c 1 1', 'I c 1 1', 'I -2xc', 'a3',
                       'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(54, 9, 'Cc', 'C c = C c 1 1', 'C c 1 1', 'C -2xc',
                       '-a1', 'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(55, 9, 'Cc', 'C c = B n 1 1', 'B n 1 1', 'B -2xbc',
                       '-a2', 'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(56, 9, 'Cc', 'C c = I b 1 1', 'I b 1 1', 'I -2xb',
                       '-a3', 'm', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(57, 10, 'P2/m', 'P 2/m = P 1 2/m 1', 'P 1 2/m 1',
                       '-P 2y', 'b', '2/m', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(58, 10, 'P2/m', 'P 2/m = P 1 1 2/m', 'P 1 1 2/m',
                       '-P 2', 'c', '2/m', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(59, 10, 'P2/m', 'P 2/m = P 2/m 1 1', 'P 2/m 1 1',
                       '-P 2x', 'a', '2/m', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(60, 11, 'P2_1/m', 'P 2_1/m = P 1 2_1/m 1',
                       'P 1 2_1/m 1', '-P 2yb', 'b', '2/m', 'monoclinic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(61, 11, 'P2_1/m', 'P 2_1/m = P 1 1 2_1/m',
                       'P 1 1 2_1/m', '-P 2c', 'c', '2/m', 'monoclinic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(62, 11, 'P2_1/m', 'P 2_1/m = P 2_1/m 1 1',
                       'P 2_1/m 1 1', '-P 2xa', 'a', '2/m', 'monoclinic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(63, 12, 'C2/m', 'C 2/m = C 1 2/m 1', 'C 1 2/m 1',
                       '-C 2y', 'b1', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(64, 12, 'C2/m', 'C 2/m = A 1 2/m 1', 'A 1 2/m 1',
                       '-A 2y', 'b2', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(65, 12, 'C2/m', 'C 2/m = I 1 2/m 1', 'I 1 2/m 1',
                       '-I 2y', 'b3', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(66, 12, 'C2/m', 'C 2/m = A 1 1 2/m', 'A 1 1 2/m',
                       '-A 2', 'c1', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(67, 12, 'C2/m', 'C 2/m = B 1 1 2/m = B 2/m',
                       'B 1 1 2/m', '-B 2', 'c2', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(68, 12, 'C2/m', 'C 2/m = I 1 1 2/m', 'I 1 1 2/m',
                       '-I 2', 'c3', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(69, 12, 'C2/m', 'C 2/m = B 2/m 1 1', 'B 2/m 1 1',
                       '-B 2x', 'a1', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(70, 12, 'C2/m', 'C 2/m = C 2/m 1 1', 'C 2/m 1 1',
                       '-C 2x', 'a2', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(71, 12, 'C2/m', 'C 2/m = I 2/m 1 1', 'I 2/m 1 1',
                       '-I 2x', 'a3', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(72, 13, 'P2/c', 'P 2/c = P 1 2/c 1', 'P 1 2/c 1',
                       '-P 2yc', 'b1', '2/m', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(73, 13, 'P2/c', 'P 2/c = P 1 2/n 1', 'P 1 2/n 1',
                       '-P 2yac', 'b2', '2/m', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(74, 13, 'P2/c', 'P 2/c = P 1 2/a 1', 'P 1 2/a 1',
                       '-P 2ya', 'b3', '2/m', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(75, 13, 'P2/c', 'P 2/c = P 1 1 2/a', 'P 1 1 2/a',
                       '-P 2a', 'c1', '2/m', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(76, 13, 'P2/c', 'P 2/c = P 1 1 2/n', 'P 1 1 2/n',
                       '-P 2ab', 'c2', '2/m', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(77, 13, 'P2/c', 'P 2/c = P 1 1 2/b = P 2/b',
                       'P 1 1 2/b', '-P 2b', 'c3', '2/m', 'monoclinic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(78, 13, 'P2/c', 'P 2/c = P 2/b 1 1', 'P 2/b 1 1',
                       '-P 2xb', 'a1', '2/m', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(79, 13, 'P2/c', 'P 2/c = P 2/n 1 1', 'P 2/n 1 1',
                       '-P 2xbc', 'a2', '2/m', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(80, 13, 'P2/c', 'P 2/c = P 2/c 1 1', 'P 2/c 1 1',
                       '-P 2xc', 'a3', '2/m', 'monoclinic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(81, 14, 'P2_1/c', 'P 2_1/c = P 1 2_1/c 1',
                       'P 1 2_1/c 1', '-P 2ybc', 'b1', '2/m', 'monoclinic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(82, 14, 'P2_1/c', 'P 2_1/c = P 1 2_1/n 1',
                       'P 1 2_1/n 1', '-P 2yn', 'b2', '2/m', 'monoclinic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(83, 14, 'P2_1/c', 'P 2_1/c = P 1 2_1/a 1',
                       'P 1 2_1/a 1', '-P 2yab', 'b3', '2/m', 'monoclinic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(84, 14, 'P2_1/c', 'P 2_1/c = P 1 1 2_1/a',
                       'P 1 1 2_1/a', '-P 2ac', 'c1', '2/m', 'monoclinic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(85, 14, 'P2_1/c', 'P 2_1/c = P 1 1 2_1/n',
                       'P 1 1 2_1/n', '-P 2n', 'c2', '2/m', 'monoclinic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(86, 14, 'P2_1/c', 'P 2_1/c = P 1 1 2_1/b = P 2_1/b',
                       'P 1 1 2_1/b', '-P 2bc', 'c3', '2/m', 'monoclinic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(87, 14, 'P2_1/c', 'P 2_1/c = P 2_1/b 1 1',
                       'P 2_1/b 1 1', '-P 2xab', 'a1', '2/m', 'monoclinic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(88, 14, 'P2_1/c', 'P 2_1/c = P 2_1/n 1 1',
                       'P 2_1/n 1 1', '-P 2xn', 'a2', '2/m', 'monoclinic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(89, 14, 'P2_1/c', 'P 2_1/c = P 2_1/c 1 1',
                       'P 2_1/c 1 1', '-P 2xac', 'a3', '2/m', 'monoclinic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(90, 15, 'C2/c', 'C 2/c = C 1 2/c 1', 'C 1 2/c 1',
                       '-C 2yc', 'b1', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(91, 15, 'C2/c', 'C 2/c = A 1 2/n 1', 'A 1 2/n 1',
                       '-A 2yac', 'b2', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(92, 15, 'C2/c', 'C 2/c = I 1 2/a 1', 'I 1 2/a 1',
                       '-I 2ya', 'b3', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(93, 15, 'C2/c', 'C 2/c = A 1 2/a 1', 'A 1 2/a 1',
                       '-A 2ya', '-b1', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(94, 15, 'C2/c', 'C 2/c = C 1 2/n 1', 'C 1 2/n 1',
                       '-C 2ybc', '-b2', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(95, 15, 'C2/c', 'C 2/c = I 1 2/c 1', 'I 1 2/c 1',
                       '-I 2yc', '-b3', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(96, 15, 'C2/c', 'C 2/c = A 1 1 2/a', 'A 1 1 2/a',
                       '-A 2a', 'c1', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(97, 15, 'C2/c', 'C 2/c = B 1 1 2/n', 'B 1 1 2/n',
                       '-B 2bc', 'c2', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(98, 15, 'C2/c', 'C 2/c = I 1 1 2/b', 'I 1 1 2/b',
                       '-I 2b', 'c3', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(99, 15, 'C2/c', 'C 2/c = B 1 1 2/b = B 2/b',
                       'B 1 1 2/b', '-B 2b', '-c1', '2/m', 'monoclinic',
                       'base'))
        self.spacegroup.append(
            SpaceGroup(100, 15, 'C2/c', 'C 2/c = A 1 1 2/n', 'A 1 1 2/n',
                       '-A 2ac', '-c2', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(101, 15, 'C2/c', 'C 2/c = I 1 1 2/a', 'I 1 1 2/a',
                       '-I 2a', '-c3', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(102, 15, 'C2/c', 'C 2/c = B 2/b 1 1', 'B 2/b 1 1',
                       '-B 2xb', 'a1', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(103, 15, 'C2/c', 'C 2/c = C 2/n 1 1', 'C 2/n 1 1',
                       '-C 2xbc', 'a2', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(104, 15, 'C2/c', 'C 2/c = I 2/c 1 1', 'I 2/c 1 1',
                       '-I 2xc', 'a3', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(105, 15, 'C2/c', 'C 2/c = C 2/c 1 1', 'C 2/c 1 1',
                       '-C 2xc', '-a1', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(106, 15, 'C2/c', 'C 2/c = B 2/n 1 1', 'B 2/n 1 1',
                       '-B 2xbc', '-a2', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(107, 15, 'C2/c', 'C 2/c = I 2/b 1 1', 'I 2/b 1 1',
                       '-I 2xb', '-a3', '2/m', 'monoclinic', 'base'))
        self.spacegroup.append(
            SpaceGroup(108, 16, 'P222', 'P 2 2 2', 'P 2 2 2', 'P 2 2', '',
                       '222', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(109, 17, 'P222_1', 'P 2 2 2_1', 'P 2 2 2_1', 'P 2c 2',
                       '', '222', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(110, 17, 'P2_122', 'P 2_1 2 2', 'P 2_1 2 2', 'P 2a 2a',
                       'cab', '222', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(111, 17, 'P22_12', 'P 2 2_1 2', 'P 2 2_1 2', 'P 2 2b',
                       'bca', '222', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(112, 18, 'P2_12_12', 'P 2_1 2_1 2', 'P 2_1 2_1 2',
                       'P 2 2ab', '', '222', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(113, 18, 'P22_12_1', 'P 2 2_1 2_1', 'P 2 2_1 2_1',
                       'P 2bc 2', 'cab', '222', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(114, 18, 'P2_122_1', 'P 2_1 2 2_1', 'P 2_1 2 2_1',
                       'P 2ac 2ac', 'bca', '222', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(115, 19, 'P2_12_12_1', 'P 2_1 2_1 2_1', 'P 2_1 2_1 2_1',
                       'P 2ac 2ab', '', '222', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(116, 20, 'C222_1', 'C 2 2 2_1', 'C 2 2 2_1', 'C 2c 2',
                       '', '222', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(117, 20, 'A2_122', 'A 2_1 2 2', 'A 2_1 2 2', 'A 2a 2a',
                       'cab', '222', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(118, 20, 'B22_12', 'B 2 2_1 2', 'B 2 2_1 2', 'B 2 2b',
                       'bca', '222', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(119, 21, 'C222', 'C 2 2 2', 'C 2 2 2', 'C 2 2', '',
                       '222', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(120, 21, 'A222', 'A 2 2 2', 'A 2 2 2', 'A 2 2', 'cab',
                       '222', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(121, 21, 'B222', 'B 2 2 2', 'B 2 2 2', 'B 2 2', 'bca',
                       '222', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(122, 22, 'F222', 'F 2 2 2', 'F 2 2 2', 'F 2 2', '',
                       '222', 'orthorhombic', 'face'))
        self.spacegroup.append(
            SpaceGroup(123, 23, 'I222', 'I 2 2 2', 'I 2 2 2', 'I 2 2', '',
                       '222', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(124, 24, 'I2_12_12_1', 'I 2_1 2_1 2_1', 'I 2_1 2_1 2_1',
                       'I 2b 2c', '', '222', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(125, 25, 'Pmm2', 'P m m 2', 'P m m 2', 'P 2 -2', '',
                       'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(126, 25, 'P2mm', 'P 2 m m', 'P 2 m m', 'P -2 2', 'cab',
                       'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(127, 25, 'Pm2m', 'P m 2 m', 'P m 2 m', 'P -2 -2', 'bca',
                       'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(128, 26, 'Pmc2_1', 'P m c 2_1', 'P m c 2_1', 'P 2c -2',
                       '', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(129, 26, 'Pcm2_1', 'P c m 2_1', 'P c m 2_1', 'P 2c -2c',
                       'ba-c', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(130, 26, 'P2_1ma', 'P 2_1 m a', 'P 2_1 m a', 'P -2a 2a',
                       'cab', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(131, 26, 'P2_1am', 'P 2_1 a m', 'P 2_1 a m', 'P -2 2a',
                       '-cba', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(132, 26, 'Pb2_1m', 'P b 2_1 m', 'P b 2_1 m', 'P -2 -2b',
                       'bca', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(133, 26, 'Pm2_1b', 'P m 2_1 b', 'P m 2_1 b', 'P -2b -2',
                       'a-cb', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(134, 27, 'Pcc2', 'P c c 2', 'P c c 2', 'P 2 -2c', '',
                       'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(135, 27, 'P2aa', 'P 2 a a', 'P 2 a a', 'P -2a 2', 'cab',
                       'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(136, 27, 'Pb2b', 'P b 2 b', 'P b 2 b', 'P -2b -2b',
                       'bca', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(137, 28, 'Pma2', 'P m a 2', 'P m a 2', 'P 2 -2a', '',
                       'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(138, 28, 'Pbm2', 'P b m 2', 'P b m 2', 'P 2 -2b',
                       'ba-c', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(139, 28, 'P2mb', 'P 2 m b', 'P 2 m b', 'P -2b 2', 'cab',
                       'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(140, 28, 'P2cm', 'P 2 c m', 'P 2 c m', 'P -2c 2',
                       '-cba', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(141, 28, 'Pc2m', 'P c 2 m', 'P c 2 m', 'P -2c -2c',
                       'bca', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(142, 28, 'Pm2a', 'P m 2 a', 'P m 2 a', 'P -2a -2a',
                       'a-cb', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(143, 29, 'Pca2_1', 'P c a 2_1', 'P c a 2_1',
                       'P 2c -2ac', '', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(144, 29, 'Pbc2_1', 'P b c 2_1', 'P b c 2_1', 'P 2c -2b',
                       'ba-c', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(145, 29, 'P2_1ab', 'P 2_1 a b', 'P 2_1 a b', 'P -2b 2a',
                       'cab', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(146, 29, 'P2_1ca', 'P 2_1 c a', 'P 2_1 c a',
                       'P -2ac 2a', '-cba', 'mm2', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(147, 29, 'Pc2_1b', 'P c 2_1 b', 'P c 2_1 b',
                       'P -2bc -2c', 'bca', 'mm2', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(148, 29, 'Pb2_1a', 'P b 2_1 a', 'P b 2_1 a',
                       'P -2a -2ab', 'a-cb', 'mm2', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(149, 30, 'Pnc2', 'P n c 2', 'P n c 2', 'P 2 -2bc', '',
                       'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(150, 30, 'Pcn2', 'P c n 2', 'P c n 2', 'P 2 -2ac',
                       'ba-c', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(151, 30, 'P2na', 'P 2 n a', 'P 2 n a', 'P -2ac 2',
                       'cab', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(152, 30, 'P2an', 'P 2 a n', 'P 2 a n', 'P -2ab 2',
                       '-cba', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(153, 30, 'Pb2n', 'P b 2 n', 'P b 2 n', 'P -2ab -2ab',
                       'bca', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(154, 30, 'Pn2b', 'P n 2 b', 'P n 2 b', 'P -2bc -2bc',
                       'a-cb', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(155, 31, 'Pmn2_1', 'P m n 2_1', 'P m n 2_1', 'P 2ac -2',
                       '', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(156, 31, 'Pnm2_1', 'P n m 2_1', 'P n m 2_1',
                       'P 2bc -2bc', 'ba-c', 'mm2', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(157, 31, 'P2_1mn', 'P 2_1 m n', 'P 2_1 m n',
                       'P -2ab 2ab', 'cab', 'mm2', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(158, 31, 'P2_1nm', 'P 2_1 n m', 'P 2_1 n m', 'P -2 2ac',
                       '-cba', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(159, 31, 'Pn2_1m', 'P n 2_1 m', 'P n 2_1 m',
                       'P -2 -2bc', 'bca', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(160, 31, 'Pm2_1n', 'P m 2_1 n', 'P m 2_1 n',
                       'P -2ab -2', 'a-cb', 'mm2', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(161, 32, 'Pba2', 'P b a 2', 'P b a 2', 'P 2 -2ab', '',
                       'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(162, 32, 'P2cb', 'P 2 c b', 'P 2 c b', 'P -2bc 2',
                       'cab', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(163, 32, 'Pc2a', 'P c 2 a', 'P c 2 a', 'P -2ac -2ac',
                       'bca', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(164, 33, 'Pna2_1', 'P n a 2_1', 'P n a 2_1', 'P 2c -2n',
                       '', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(165, 33, 'Pbn2_1', 'P b n 2_1', 'P b n 2_1',
                       'P 2c -2ab', 'ba-c', 'mm2', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(166, 33, 'P2_1nb', 'P 2_1 n b', 'P 2_1 n b',
                       'P -2bc 2a', 'cab', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(167, 33, 'P2_1cn', 'P 2_1 c n', 'P 2_1 c n', 'P -2n 2a',
                       '-cba', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(168, 33, 'Pc2_1n', 'P c 2_1 n', 'P c 2_1 n',
                       'P -2n -2ac', 'bca', 'mm2', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(169, 33, 'Pn2_1a', 'P n 2_1 a', 'P n 2_1 a',
                       'P -2ac -2n', 'a-cb', 'mm2', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(170, 34, 'Pnn2', 'P n n 2', 'P n n 2', 'P 2 -2n', '',
                       'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(171, 34, 'P2nn', 'P 2 n n', 'P 2 n n', 'P -2n 2', 'cab',
                       'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(172, 34, 'Pn2n', 'P n 2 n', 'P n 2 n', 'P -2n -2n',
                       'bca', 'mm2', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(173, 35, 'Cmm2', 'C m m 2', 'C m m 2', 'C 2 -2', '',
                       'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(174, 35, 'A2mm', 'A 2 m m', 'A 2 m m', 'A -2 2', 'cab',
                       'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(175, 35, 'Bm2m', 'B m 2 m', 'B m 2 m', 'B -2 -2', 'bca',
                       'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(176, 36, 'Cmc2_1', 'C m c 2_1', 'C m c 2_1', 'C 2c -2',
                       '', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(177, 36, 'Ccm2_1', 'C c m 2_1', 'C c m 2_1', 'C 2c -2c',
                       'ba-c', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(178, 36, 'A2_1ma', 'A 2_1 m a', 'A 2_1 m a', 'A -2a 2a',
                       'cab', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(179, 36, 'A2_1am', 'A 2_1 a m', 'A 2_1 a m', 'A -2 2a',
                       '-cba', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(180, 36, 'Bb2_1m', 'B b 2_1 m', 'B b 2_1 m', 'B -2 -2b',
                       'bca', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(181, 36, 'Bm2_1b', 'B m 2_1 b', 'B m 2_1 b', 'B -2b -2',
                       'a-cb', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(182, 37, 'Ccc2', 'C c c 2', 'C c c 2', 'C 2 -2c', '',
                       'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(183, 37, 'A2aa', 'A 2 a a', 'A 2 a a', 'A -2a 2', 'cab',
                       'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(184, 37, 'Bb2b', 'B b 2 b', 'B b 2 b', 'B -2b -2b',
                       'bca', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(185, 38, 'Amm2', 'A m m 2', 'A m m 2', 'A 2 -2', '',
                       'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(186, 38, 'Bmm2', 'B m m 2', 'B m m 2', 'B 2 -2', 'ba-c',
                       'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(187, 38, 'B2mm', 'B 2 m m', 'B 2 m m', 'B -2 2', 'cab',
                       'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(188, 38, 'C2mm', 'C 2 m m', 'C 2 m m', 'C -2 2', '-cba',
                       'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(189, 38, 'Cm2m', 'C m 2 m', 'C m 2 m', 'C -2 -2', 'bca',
                       'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(190, 38, 'Am2m', 'A m 2 m', 'A m 2 m', 'A -2 -2',
                       'a-cb', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(191, 39, 'Aem2', 'A e m 2', 'A e m 2', 'A 2 -2c', '',
                       'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(192, 39, 'Bme2', 'B m e 2', 'B m e 2', 'B 2 -2c',
                       'ba-c', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(193, 39, 'B2em', 'B 2 e m', 'B 2 e m', 'B -2c 2', 'cab',
                       'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(194, 39, 'C2me', 'C 2 m e', 'C 2 m e', 'C -2b 2',
                       '-cba', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(195, 39, 'Cm2e', 'C m 2 e', 'C m 2 e', 'C -2b -2b',
                       'bca', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(196, 39, 'Ae2m', 'A e 2 m', 'A e 2 m', 'A -2c -2c',
                       'a-cb', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(197, 40, 'Ama2', 'A m a 2', 'A m a 2', 'A 2 -2a', '',
                       'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(198, 40, 'Bbm2', 'B b m 2', 'B b m 2', 'B 2 -2b',
                       'ba-c', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(199, 40, 'B2mb', 'B 2 m b', 'B 2 m b', 'B -2b 2', 'cab',
                       'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(200, 40, 'C2cm', 'C 2 c m', 'C 2 c m', 'C -2c 2',
                       '-cba', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(201, 40, 'Cc2m', 'C c 2 m', 'C c 2 m', 'C -2c -2c',
                       'bca', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(202, 40, 'Am2a', 'A m 2 a', 'A m 2 a', 'A -2a -2a',
                       'a-cb', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(203, 41, 'Aea2', 'A e a 2', 'A e a 2', 'A 2 -2ac', '',
                       'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(204, 41, 'Bbe2', 'B b e 2', 'B b e 2', 'B 2 -2bc',
                       'ba-c', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(205, 41, 'B2eb', 'B 2 e b', 'B 2 e b', 'B -2bc 2',
                       'cab', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(206, 41, 'C2ce', 'C 2 c e', 'C 2 c e', 'C -2bc 2',
                       '-cba', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(207, 41, 'Cc2e', 'C c 2 e', 'C c 2 e', 'C -2bc -2bc',
                       'bca', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(208, 41, 'Ae2a', 'A e 2 a', 'A e 2 a', 'A -2ac -2ac',
                       'a-cb', 'mm2', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(209, 42, 'Fmm2', 'F m m 2', 'F m m 2', 'F 2 -2', '',
                       'mm2', 'orthorhombic', 'face'))
        self.spacegroup.append(
            SpaceGroup(210, 42, 'F2mm', 'F 2 m m', 'F 2 m m', 'F -2 2', 'cab',
                       'mm2', 'orthorhombic', 'face'))
        self.spacegroup.append(
            SpaceGroup(211, 42, 'Fm2m', 'F m 2 m', 'F m 2 m', 'F -2 -2', 'bca',
                       'mm2', 'orthorhombic', 'face'))
        self.spacegroup.append(
            SpaceGroup(212, 43, 'Fdd2', 'F d d 2', 'F d d 2', 'F 2 -2d', '',
                       'mm2', 'orthorhombic', 'face'))
        self.spacegroup.append(
            SpaceGroup(213, 43, 'F2dd', 'F 2 d d', 'F 2 d d', 'F -2d 2', 'cab',
                       'mm2', 'orthorhombic', 'face'))
        self.spacegroup.append(
            SpaceGroup(214, 43, 'Fd2d', 'F d 2 d', 'F d 2 d', 'F -2d -2d',
                       'bca', 'mm2', 'orthorhombic', 'face'))
        self.spacegroup.append(
            SpaceGroup(215, 44, 'Imm2', 'I m m 2', 'I m m 2', 'I 2 -2', '',
                       'mm2', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(216, 44, 'I2mm', 'I 2 m m', 'I 2 m m', 'I -2 2', 'cab',
                       'mm2', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(217, 44, 'Im2m', 'I m 2 m', 'I m 2 m', 'I -2 -2', 'bca',
                       'mm2', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(218, 45, 'Iba2', 'I b a 2', 'I b a 2', 'I 2 -2c', '',
                       'mm2', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(219, 45, 'I2cb', 'I 2 c b', 'I 2 c b', 'I -2a 2', 'cab',
                       'mm2', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(220, 45, 'Ic2a', 'I c 2 a', 'I c 2 a', 'I -2b -2b',
                       'bca', 'mm2', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(221, 46, 'Ima2', 'I m a 2', 'I m a 2', 'I 2 -2a', '',
                       'mm2', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(222, 46, 'Ibm2', 'I b m 2', 'I b m 2', 'I 2 -2b',
                       'ba-c', 'mm2', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(223, 46, 'I2mb', 'I 2 m b', 'I 2 m b', 'I -2b 2', 'cab',
                       'mm2', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(224, 46, 'I2cm', 'I 2 c m', 'I 2 c m', 'I -2c 2',
                       '-cba', 'mm2', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(225, 46, 'Ic2m', 'I c 2 m', 'I c 2 m', 'I -2c -2c',
                       'bca', 'mm2', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(226, 46, 'Im2a', 'I m 2 a', 'I m 2 a', 'I -2a -2a',
                       'a-cb', 'mm2', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(227, 47, 'Pmmm', 'P m m m', 'P 2/m 2/m 2/m', '-P 2 2',
                       '', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(228, 48, 'Pnnn', 'P n n n', 'P 2/n 2/n 2/n',
                       'P 2 2 -1n', '1', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(229, 48, 'Pnnn', 'P n n n', 'P 2/n 2/n 2/n',
                       '-P 2ab 2bc', '2', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(230, 49, 'Pccm', 'P c c m', 'P 2/c 2/c 2/m', '-P 2 2c',
                       '', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(231, 49, 'Pmaa', 'P m a a', 'P 2/m 2/a 2/a', '-P 2a 2',
                       'cab', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(232, 49, 'Pbmb', 'P b m b', 'P 2/b 2/m 2/b', '-P 2b 2b',
                       'bca', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(233, 50, 'Pban', 'P b a n', 'P 2/b 2/a 2/n',
                       'P 2 2 -1ab', '1', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(234, 50, 'Pban', 'P b a n', 'P 2/b 2/a 2/n',
                       '-P 2ab 2b', '2', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(235, 50, 'Pncb', 'P n c b', 'P 2/n 2/c 2/b',
                       'P 2 2 -1bc', '1cab', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(236, 50, 'Pncb', 'P n c b', 'P 2/n 2/c 2/b',
                       '-P 2b 2bc', '2cab', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(237, 50, 'Pcna', 'P c n a', 'P 2/c 2/n 2/a',
                       'P 2 2 -1ac', '1bca', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(238, 50, 'Pcna', 'P c n a', 'P 2/c 2/n 2/a', '-P 2a 2c',
                       '2bca', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(239, 51, 'Pmma', 'P m m a', 'P 2_1/m 2/m 2/a',
                       '-P 2a 2a', '', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(240, 51, 'Pmmb', 'P m m b', 'P 2/m 2_1/m 2/b',
                       '-P 2b 2', 'ba-c', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(241, 51, 'Pbmm', 'P b m m', 'P 2/b 2_1/m 2/m',
                       '-P 2 2b', 'cab', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(242, 51, 'Pcmm', 'P c m m', 'P 2/c 2/m 2_1/m',
                       '-P 2c 2c', '-cba', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(243, 51, 'Pmcm', 'P m c m', 'P 2/m 2/c 2_1/m',
                       '-P 2c 2', 'bca', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(244, 51, 'Pmam', 'P m a m', 'P 2_1/m 2/a 2/m',
                       '-P 2 2a', 'a-cb', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(245, 52, 'Pnna', 'P n n a', 'P 2/n 2_1/n 2/a',
                       '-P 2a 2bc', '', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(246, 52, 'Pnnb', 'P n n b', 'P 2_1/n 2/n 2/b',
                       '-P 2b 2n', 'ba-c', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(247, 52, 'Pbnn', 'P b n n', 'P 2/b 2/n 2_1/n',
                       '-P 2n 2b', 'cab', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(248, 52, 'Pcnn', 'P c n n', 'P 2/c 2_1/n 2/n',
                       '-P 2ab 2c', '-cba', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(249, 52, 'Pncn', 'P n c n', 'P 2_1/n 2/c 2/n',
                       '-P 2ab 2n', 'bca', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(250, 52, 'Pnan', 'P n a n', 'P 2/n 2/a 2_1/n',
                       '-P 2n 2bc', 'a-cb', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(251, 53, 'Pmna', 'P m n a', 'P 2/m 2/n 2_1/a',
                       '-P 2ac 2', '', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(252, 53, 'Pnmb', 'P n m b', 'P 2/n 2/m 2_1/b',
                       '-P 2bc 2bc', 'ba-c', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(253, 53, 'Pbmn', 'P b m n', 'P 2_1/b 2/m 2/n',
                       '-P 2ab 2ab', 'cab', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(254, 53, 'Pcnm', 'P c n m', 'P 2_1/c 2/n 2/m',
                       '-P 2 2ac', '-cba', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(255, 53, 'Pncm', 'P n c m', 'P 2/n 2_1/c 2/m',
                       '-P 2 2bc', 'bca', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(256, 53, 'Pman', 'P m a n', 'P 2/m 2_1/a 2/n',
                       '-P 2ab 2', 'a-cb', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(257, 54, 'Pcca', 'P c c a', 'P 2_1/c 2/c 2/a',
                       '-P 2a 2ac', '', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(258, 54, 'Pccb', 'P c c b', 'P 2/c 2_1/c 2/b',
                       '-P 2b 2c', 'ba-c', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(259, 54, 'Pbaa', 'P b a a', 'P 2/b 2_1/a 2/a',
                       '-P 2a 2b', 'cab', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(260, 54, 'Pcaa', 'P c a a', 'P 2/c 2/a 2_1/a',
                       '-P 2ac 2c', '-cba', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(261, 54, 'Pbcb', 'P b c b', 'P 2/b 2/c 2_1/b',
                       '-P 2bc 2b', 'bca', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(262, 54, 'Pbab', 'P b a b', 'P 2_1/b 2/a 2/b',
                       '-P 2b 2ab', 'a-cb', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(263, 55, 'Pbam', 'P b a m', 'P 2_1/b 2_1/a 2/m',
                       '-P 2 2ab', '', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(264, 55, 'Pmcb', 'P m c b', 'P 2/m 2_1/c 2_1/b',
                       '-P 2bc 2', 'cab', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(265, 55, 'Pcma', 'P c m a', 'P 2_1/c 2/m 2_1/a',
                       '-P 2ac 2ac', 'bca', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(266, 56, 'Pccn', 'P c c n', 'P 2_1/c 2_1/c 2/n',
                       '-P 2ab 2ac', '', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(267, 56, 'Pnaa', 'P n a a', 'P 2/n 2_1/a 2_1/a',
                       '-P 2ac 2bc', 'cab', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(268, 56, 'Pbnb', 'P b n b', 'P 2_1/b 2/n 2_1/b',
                       '-P 2bc 2ab', 'bca', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(269, 57, 'Pbcm', 'P b c m', 'P 2/b 2_1/c 2_1/m',
                       '-P 2c 2b', '', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(270, 57, 'Pcam', 'P c a m', 'P 2_1/c 2/a 2_1/m',
                       '-P 2c 2ac', 'ba-c', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(271, 57, 'Pmca', 'P m c a', 'P 2_1/m 2/c 2_1/a',
                       '-P 2ac 2a', 'cab', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(272, 57, 'Pmab', 'P m a b', 'P 2_1/m 2_1/a 2/b',
                       '-P 2b 2a', '-cba', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(273, 57, 'Pbma', 'P b m a', 'P 2_1/b 2_1/m 2/a',
                       '-P 2a 2ab', 'bca', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(274, 57, 'Pcmb', 'P c m b', 'P 2/c 2_1/m 2_1/b',
                       '-P 2bc 2c', 'a-cb', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(275, 58, 'Pnnm', 'P n n m', 'P 2_1/n 2_1/n 2/m',
                       '-P 2 2n', '', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(276, 58, 'Pmnn', 'P m n n', 'P 2/m 2_1/n 2_1/n',
                       '-P 2n 2', 'cab', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(277, 58, 'Pnmn', 'P n m n', 'P 2_1/n 2/m 2_1/n',
                       '-P 2n 2n', 'bca', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(278, 59, 'Pmmn', 'P m m n', 'P 2_1/m 2_1/m 2/n',
                       'P 2 2ab -1ab', '1', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(279, 59, 'Pmmn', 'P m m n', 'P 2_1/m 2_1/m 2/n',
                       '-P 2ab 2a', '2', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(280, 59, 'Pnmm', 'P n m m', 'P 2/n 2_1/m 2_1/m',
                       'P 2bc 2 -1bc', '1cab', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(281, 59, 'Pnmm', 'P n m m', 'P 2/n 2_1/m 2_1/m',
                       '-P 2c 2bc', '2cab', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(282, 59, 'Pmnm', 'P m n m', 'P 2_1/m 2/n 2_1/m',
                       'P 2ac 2ac -1ac', '1bca', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(283, 59, 'Pmnm', 'P m n m', 'P 2_1/m 2/n 2_1/m',
                       '-P 2c 2a', '2bca', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(284, 60, 'Pbcn', 'P b c n', 'P 2_1/b 2/c 2_1/n',
                       '-P 2n 2ab', '', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(285, 60, 'Pcan', 'P c a n', 'P 2/c 2_1/a 2_1/n',
                       '-P 2n 2c', 'ba-c', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(286, 60, 'Pnca', 'P n c a', 'P 2_1/n 2_1/c 2/a',
                       '-P 2a 2n', 'cab', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(287, 60, 'Pnab', 'P n a b', 'P 2_1/n 2/a 2_1/b',
                       '-P 2bc 2n', '-cba', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(288, 60, 'Pbna', 'P b n a', 'P 2/b 2_1/n 2_1/a',
                       '-P 2ac 2b', 'bca', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(289, 60, 'Pcnb', 'P c n b', 'P 2_1/c 2_1/n 2/b',
                       '-P 2b 2ac', 'a-cb', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(290, 61, 'Pbca', 'P b c a', 'P 2_1/b 2_1/c 2_1/a',
                       '-P 2ac 2ab', '', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(291, 61, 'Pcab', 'P c a b', 'P 2_1/c 2_1/a 2_1/b',
                       '-P 2bc 2ac', 'ba-c', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(292, 62, 'Pnma', 'P n m a', 'P 2_1/n 2_1/m 2_1/a',
                       '-P 2ac 2n', '', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(293, 62, 'Pmnb', 'P m n b', 'P 2_1/m 2_1/n 2_1/b',
                       '-P 2bc 2a', 'ba-c', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(294, 62, 'Pbnm', 'P b n m', 'P 2_1/b 2_1/n 2_1/m',
                       '-P 2c 2ab', 'cab', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(295, 62, 'Pcmn', 'P c m n', 'P 2_1/c 2_1/m 2_1/n',
                       '-P 2n 2ac', '-cba', 'mmm', 'orthorhombic',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(296, 62, 'Pmcn', 'P m c n', 'P 2_1/m 2_1/c 2_1/n',
                       '-P 2n 2a', 'bca', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(297, 62, 'Pnam', 'P n a m', 'P 2_1/n 2_1/a 2_1/m',
                       '-P 2c 2n', 'a-cb', 'mmm', 'orthorhombic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(298, 63, 'Cmcm', 'C m c m', 'C 2/m 2/c 2_1/m',
                       '-C 2c 2', '', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(299, 63, 'Ccmm', 'C c m m', 'C 2/c 2/m 2_1/m',
                       '-C 2c 2c', 'ba-c', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(300, 63, 'Amma', 'A m m a', 'A 2_1/m 2/m 2/a',
                       '-A 2a 2a', 'cab', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(301, 63, 'Amam', 'A m a m', 'A 2_1/m 2/a 2/m',
                       '-A 2 2a', '-cba', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(302, 63, 'Bbmm', 'B b m m', 'B 2/b 2_1/m 2/m',
                       '-B 2 2b', 'bca', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(303, 63, 'Bmmb', 'B m m b', 'B 2/m 2_1/m 2/b',
                       '-B 2b 2', 'a-cb', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(304, 64, 'Cmce', 'C m c e', 'C 2/m 2/c 2_1/e',
                       '-C 2bc 2', '', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(305, 64, 'Ccme', 'C c m e', 'C 2/c 2/m 2_1/e',
                       '-C 2bc 2bc', 'ba-c', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(306, 64, 'Aema', 'A e m a', 'A 2_1/e 2/m 2/a',
                       '-A 2ac 2ac', 'cab', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(307, 64, 'Aeam', 'A e a m', 'A 2_1/e 2/a 2/m',
                       '-A 2 2ac', '-cba', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(308, 64, 'Bbem', 'B b e m', 'B 2/b 2_1/e 2/m',
                       '-B 2 2bc', 'bca', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(309, 64, 'Bmeb', 'B m e b', 'B 2/m 2_1/e 2/b',
                       '-B 2bc 2', 'a-cb', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(310, 65, 'Cmmm', 'C m m m', 'C 2/m 2/m 2/m', '-C 2 2',
                       '', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(311, 65, 'Ammm', 'A m m m', 'A 2/m 2/m 2/m', '-A 2 2',
                       'cab', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(312, 65, 'Bmmm', 'B m m m', 'B 2/m 2/m 2/m', '-B 2 2',
                       'bca', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(313, 66, 'Cccm', 'C c c m', 'C 2/c 2/c 2/m', '-C 2 2c',
                       '', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(314, 66, 'Amaa', 'A m a a', 'A 2/m 2/a 2/a', '-A 2a 2',
                       'cab', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(315, 66, 'Bbmb', 'B b m b', 'B 2/b 2/m 2/b', '-B 2b 2b',
                       'bca', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(316, 67, 'Cmme', 'C m m e', 'C 2/m 2/m 2/e', '-C 2b 2',
                       '', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(317, 67, 'Cmme', 'C m m e', 'C 2/m 2/m 2/e', '-C 2b 2b',
                       'ba-c', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(318, 67, 'Aemm', 'A e m m', 'A 2/e 2/m 2/m', '-A 2c 2c',
                       'cab', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(319, 67, 'Aemm', 'A e m m', 'A 2/e 2/m 2/m', '-A 2 2c',
                       '-cba', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(320, 67, 'Bmem', 'B m e m', 'B 2/m 2/e 2/m', '-B 2 2c',
                       'bca', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(321, 67, 'Bmem', 'B m e m', 'B 2/m 2/e 2/m', '-B 2c 2',
                       'a-cb', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(322, 68, 'Ccce', 'C c c e', 'C 2/c 2/c 2/e',
                       'C 2 2 -1bc', '1', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(323, 68, 'Ccce', 'C c c e', 'C 2/c 2/c 2/e',
                       '-C 2b 2bc', '2', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(324, 68, 'Ccce', 'C c c e', 'C 2/c 2/c 2/e',
                       'C 2 2 -1bc', '1ba-c', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(325, 68, 'Ccce', 'C c c e', 'C 2/c 2/c 2/e', '-C 2b 2c',
                       '2ba-c', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(326, 68, 'Aeaa', 'A e a a', 'A 2/e 2/a 2/a',
                       'A 2 2 -1ac', '1cab', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(327, 68, 'Aeaa', 'A e a a', 'A 2/e 2/a 2/a', '-A 2a 2c',
                       '2cab', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(328, 68, 'Aeaa', 'A e a a', 'A 2/e 2/a 2/a',
                       'A 2 2 -1ac', '1-cba', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(329, 68, 'Aeaa', 'A e a a', 'A 2/e 2/a 2/a',
                       '-A 2ac 2c', '2-cba', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(330, 68, 'Bbeb', 'B b e b', 'B 2/b 2/e 2/b',
                       'B 2 2 -1bc', '1bca', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(331, 68, 'Bbcb', 'B b c b', 'B 2/b 2/e 2/b',
                       '-B 2bc 2b', '2bca', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(332, 68, 'Bbeb', 'B b e b', 'B 2/b 2/e 2/b',
                       'B 2 2 -1bc', '1a-cb', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(333, 68, 'Bbeb', 'B b e b', 'B 2/b 2/e 2/b',
                       '-B 2b 2bc', '2a-cb', 'mmm', 'orthorhombic', 'base'))
        self.spacegroup.append(
            SpaceGroup(334, 69, 'Fmmm', 'F m m m', 'F 2/m 2/m 2/m', '-F 2 2',
                       '', 'mmm', 'orthorhombic', 'face'))
        self.spacegroup.append(
            SpaceGroup(335, 70, 'Fddd', 'F d d d', 'F 2/d 2/d 2/d',
                       'F 2 2 -1d', '1', 'mmm', 'orthorhombic', 'face'))
        self.spacegroup.append(
            SpaceGroup(336, 70, 'Fddd', 'F d d d', 'F 2/d 2/d 2/d',
                       '-F 2uv 2vw', '2', 'mmm', 'orthorhombic', 'face'))
        self.spacegroup.append(
            SpaceGroup(337, 71, 'Immm', 'I m m m', 'I 2/m 2/m 2/m', '-I 2 2',
                       '', 'mmm', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(338, 72, 'Ibam', 'I b a m', 'I 2/b 2/a 2/m', '-I 2 2c',
                       '', 'mmm', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(339, 72, 'Imcb', 'I m c b', 'I 2/m 2/c 2/b', '-I 2a 2',
                       'cab', 'mmm', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(340, 72, 'Icma', 'I c m a', 'I 2/c 2/m 2/a', '-I 2b 2b',
                       'bca', 'mmm', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(341, 73, 'Ibca', 'I b c a', 'I 2/b 2/c 2/a', '-I 2b 2c',
                       '', 'mmm', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(342, 73, 'Icab', 'I c a b', 'I 2/c 2/a 2/b', '-I 2a 2b',
                       'ba-c', 'mmm', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(343, 74, 'Imma', 'I m m a', 'I 2/m 2/m 2/a', '-I 2b 2',
                       '', 'mmm', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(344, 74, 'Immb', 'I m m b', 'I 2/m 2/m 2/b', '-I 2a 2a',
                       'ba-c', 'mmm', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(345, 74, 'Ibmm', 'I b m m', 'I 2/b 2/m 2/m', '-I 2c 2c',
                       'cab', 'mmm', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(346, 74, 'Icmm', 'I c m m', 'I 2/c 2/m 2/m', '-I 2 2b',
                       '-cba', 'mmm', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(347, 74, 'Imcm', 'I m c m', 'I 2/m 2/c 2/m', '-I 2 2a',
                       'bca', 'mmm', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(348, 74, 'Imam', 'I m a m', 'I 2/m 2/a 2/m', '-I 2c 2',
                       'a-cb', 'mmm', 'orthorhombic', 'body'))
        self.spacegroup.append(
            SpaceGroup(349, 75, 'P4', 'P 4', 'P 4', 'P 4', '', '4',
                       'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(350, 76, 'P4_1', 'P 4_1', 'P 4_1', 'P 4w', '', '4',
                       'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(351, 77, 'P4_2', 'P 4_2', 'P 4_2', 'P 4c', '', '4',
                       'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(352, 78, 'P4_3', 'P 4_3', 'P 4_3', 'P 4cw', '', '4',
                       'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(353, 79, 'I4', 'I 4', 'I 4', 'I 4', '', '4',
                       'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(354, 80, 'I4_1', 'I 4_1', 'I 4_1', 'I 4bw', '', '4',
                       'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(355, 81, 'P-4', 'P -4', 'P -4', 'P -4', '', '-4',
                       'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(356, 82, 'I-4', 'I -4', 'I -4', 'I -4', '', '-4',
                       'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(357, 83, 'P4/m', 'P 4/m', 'P 4/m', '-P 4', '', '4/m',
                       'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(358, 84, 'P4_2/m', 'P 4_2/m', 'P 4_2/m', '-P 4c', '',
                       '4/m', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(359, 85, 'P4/n', 'P 4/n', 'P 4/n', 'P 4ab -1ab', '1',
                       '4/m', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(360, 85, 'P4/n', 'P 4/n', 'P 4/n', '-P 4a', '2', '4/m',
                       'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(361, 86, 'P4_2/n', 'P 4_2/n', 'P 4_2/n', 'P 4n -1n',
                       '1', '4/m', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(362, 86, 'P4_2/n', 'P 4_2/n', 'P 4_2/n', '-P 4bc', '2',
                       '4/m', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(363, 87, 'I4/m', 'I 4/m', 'I 4/m', '-I 4', '', '4/m',
                       'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(364, 88, 'I4_1/a', 'I 4_1/a', 'I 4_1/a', 'I 4bw -1bw',
                       '1', '4/m', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(365, 88, 'I4_1/a', 'I 4_1/a', 'I 4_1/a', '-I 4ad', '2',
                       '4/m', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(366, 89, 'P422', 'P 4 2 2', 'P 4 2 2', 'P 4 2', '',
                       '422', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(367, 90, 'P42_12', 'P 4 2_1 2', 'P 4 2_1 2',
                       'P 4ab 2ab', '', '422', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(368, 91, 'P4_122', 'P 4_1 2 2', 'P 4_1 2 2', 'P 4w 2c',
                       '', '422', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(369, 92, 'P4_12_12', 'P 4_1 2_1 2', 'P 4_1 2_1 2',
                       'P 4abw 2nw', '', '422', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(370, 93, 'P4_222', 'P 4_2 2 2', 'P 4_2 2 2', 'P 4c 2',
                       '', '422', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(371, 94, 'P4_22_12', 'P 4_2 2_1 2', 'P 4_2 2_1 2',
                       'P 4n 2n', '', '422', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(372, 95, 'P4_322', 'P 4_3 2 2', 'P 4_3 2 2', 'P 4cw 2c',
                       '', '422', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(373, 96, 'P4_32_12', 'P 4_3 2_1 2', 'P 4_3 2_1 2',
                       'P 4nw 2abw', '', '422', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(374, 97, 'I422', 'I 4 2 2', 'I 4 2 2', 'I 4 2', '',
                       '422', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(375, 98, 'I4_122', 'I 4_1 2 2', 'I 4_1 2 2',
                       'I 4bw 2bw', '', '422', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(376, 99, 'P4mm', 'P 4 m m', 'P 4 m m', 'P 4 -2', '',
                       '4mm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(377, 100, 'P4bm', 'P 4 b m', 'P 4 b m', 'P 4 -2ab', '',
                       '4mm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(378, 101, 'P4_2cm', 'P 4_2 c m', 'P 4_2 c m',
                       'P 4c -2c', '', '4mm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(379, 102, 'P4_2nm', 'P 4_2 n m', 'P 4_2 n m',
                       'P 4n -2n', '', '4mm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(380, 103, 'P4cc', 'P 4 c c', 'P 4 c c', 'P 4 -2c', '',
                       '4mm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(381, 104, 'P4nc', 'P 4 n c', 'P 4 n c', 'P 4 -2n', '',
                       '4mm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(382, 105, 'P4_2mc', 'P 4_2 m c', 'P 4_2 m c', 'P 4c -2',
                       '', '4mm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(383, 106, 'P4_2bc', 'P 4_2 b c', 'P 4_2 b c',
                       'P 4c -2ab', '', '4mm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(384, 107, 'I4mm', 'I 4 m m', 'I 4 m m', 'I 4 -2', '',
                       '4mm', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(385, 108, 'I4cm', 'I 4 c m', 'I 4 c m', 'I 4 -2c', '',
                       '4mm', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(386, 109, 'I4_1md', 'I 4_1 m d', 'I 4_1 m d',
                       'I 4bw -2', '', '4mm', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(387, 110, 'I4_1cd', 'I 4_1 c d', 'I 4_1 c d',
                       'I 4bw -2c', '', '4mm', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(388, 111, 'P-42m', 'P -4 2 m', 'P -4 2 m', 'P -4 2', '',
                       '-42m', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(389, 112, 'P-42c', 'P -4 2 c', 'P -4 2 c', 'P -4 2c',
                       '', '-42m', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(390, 113, 'P-42_1m', 'P -4 2_1 m', 'P -4 2_1 m',
                       'P -4 2ab', '', '-42m', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(391, 114, 'P-42_1c', 'P -4 2_1 c', 'P -4 2_1 c',
                       'P -4 2n', '', '-42m', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(392, 115, 'P-4m2', 'P -4 m 2', 'P -4 m 2', 'P -4 -2',
                       '', '-42m', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(393, 116, 'P-4c2', 'P -4 c 2', 'P -4 c 2', 'P -4 -2c',
                       '', '-42m', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(394, 117, 'P-4b2', 'P -4 b 2', 'P -4 b 2', 'P -4 -2ab',
                       '', '-42m', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(395, 118, 'P-4n2', 'P -4 n 2', 'P -4 n 2', 'P -4 -2n',
                       '', '-42m', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(396, 119, 'I-4m2', 'I -4 m 2', 'I -4 m 2', 'I -4 -2',
                       '', '-42m', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(397, 120, 'I-4c2', 'I -4 c 2', 'I -4 c 2', 'I -4 -2c',
                       '', '-42m', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(398, 121, 'I-42m', 'I -4 2 m', 'I -4 2 m', 'I -4 2', '',
                       '-42m', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(399, 122, 'I-42d', 'I -4 2 d', 'I -4 2 d', 'I -4 2bw',
                       '', '-42m', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(400, 123, 'P4/mmm', 'P 4/m m m', 'P 4/m 2/m 2/m',
                       '-P 4 2', '', '4/mmm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(401, 124, 'P4/mcc', 'P 4/m c c', 'P 4/m 2/c 2/c',
                       '-P 4 2c', '', '4/mmm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(402, 125, 'P4/nbm', 'P 4/n b m', 'P 4/n 2/b 2/m',
                       'P 4 2 -1ab', '1', '4/mmm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(403, 125, 'P4/nbm', 'P 4/n b m', 'P 4/n 2/b 2/m',
                       '-P 4a 2b', '2', '4/mmm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(404, 126, 'P4/nnc', 'P 4/n n c', 'P 4/n 2/n 2/c',
                       'P 4 2 -1n', '1', '4/mmm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(405, 126, 'P4/nnc', 'P 4/n n c', 'P 4/n 2/n 2/c',
                       '-P 4a 2bc', '2', '4/mmm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(406, 127, 'P4/mbm', 'P 4/m b m', 'P 4/m 2_1/b m',
                       '-P 4 2ab', '', '4/mmm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(407, 128, 'P4/mnc', 'P 4/m n c', 'P 4/m 2_1/n c',
                       '-P 4 2n', '', '4/mmm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(408, 129, 'P4/nmm', 'P 4/n m m', 'P 4/n 2_1/m m',
                       'P 4ab 2ab -1ab', '1', '4/mmm', 'tetragonal',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(409, 129, 'P4/nmm', 'P 4/n m m', 'P 4/n 2_1/m m',
                       '-P 4a 2a', '2', '4/mmm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(410, 130, 'P4/ncc', 'P 4/n c c', 'P 4/n 2_1/c c',
                       'P 4ab 2n -1ab', '1', '4/mmm', 'tetragonal',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(411, 130, 'P4/ncc', 'P 4/n c c', 'P 4/n 2_1/c c',
                       '-P 4a 2ac', '2', '4/mmm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(412, 131, 'P4_2/mmc', 'P 4_2/m m c', 'P 4_2/m 2/m 2/c',
                       '-P 4c 2', '', '4/mmm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(413, 132, 'P4_2/mcm', 'P 4_2/m c m', 'P 4_2/m 2/c 2/m',
                       '-P 4c 2c', '', '4/mmm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(414, 133, 'P4_2/nbc', 'P 4_2/n b c', 'P 4_2/n 2/b 2/c',
                       'P 4n 2c -1n', '1', '4/mmm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(415, 133, 'P4_2/nbc', 'P 4_2/n b c', 'P 4_2/n 2/b 2/c',
                       '-P 4ac 2b', '2', '4/mmm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(416, 134, 'P4_2/nnm', 'P 4_2/n n m', 'P 4_2/n 2/n 2/m',
                       'P 4n 2 -1n', '1', '4/mmm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(417, 134, 'P4_2/nnm', 'P 4_2/n n m', 'P 4_2/n 2/n 2/m',
                       '-P 4ac 2bc', '2', '4/mmm', 'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(418, 135, 'P4_2/mbc', 'P 4_2/m b c',
                       'P 4_2/m 2_1/b 2/c', '-P 4c 2ab', '', '4/mmm',
                       'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(419, 136, 'P4_2/mnm', 'P 4_2/m n m',
                       'P 4_2/m 2_1/n 2/m', '-P 4n 2n', '', '4/mmm',
                       'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(420, 137, 'P4_2/nmc', 'P 4_2/n m c',
                       'P 4_2/n 2_1/m 2/c', 'P 4n 2n -1n', '1', '4/mmm',
                       'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(421, 137, 'P4_2/nmc', 'P 4_2/n m c',
                       'P 4_2/n 2_1/m 2/c', '-P 4ac 2a', '2', '4/mmm',
                       'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(422, 138, 'P4_2/ncm', 'P 4_2/n c m',
                       'P 4_2/n 2_1/c 2/m', 'P 4n 2ab -1n', '1', '4/mmm',
                       'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(423, 138, 'P4_2/ncm', 'P 4_2/n c m',
                       'P 4_2/n 2_1/c 2/m', '-P 4ac 2ac', '2', '4/mmm',
                       'tetragonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(424, 139, 'I4/mmm', 'I 4/m m m', 'I 4/m 2/m 2/m',
                       '-I 4 2', '', '4/mmm', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(425, 140, 'I4/mcm', 'I 4/m c m', 'I 4/m 2/c 2/m',
                       '-I 4 2c', '', '4/mmm', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(426, 141, 'I4_1/amd', 'I 4_1/a m d', 'I 4_1/a 2/m 2/d',
                       'I 4bw 2bw -1bw', '1', '4/mmm', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(427, 141, 'I4_1/amd', 'I 4_1/a m d', 'I 4_1/a 2/m 2/d',
                       '-I 4bd 2', '2', '4/mmm', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(428, 142, 'I4_1/acd', 'I 4_1/a c d', 'I 4_1/a 2/c 2/d',
                       'I 4bw 2aw -1bw', '1', '4/mmm', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(429, 142, 'I4_1/acd', 'I 4_1/a c d', 'I 4_1/a 2/c 2/d',
                       '-I 4bd 2c', '2', '4/mmm', 'tetragonal', 'body'))
        self.spacegroup.append(
            SpaceGroup(430, 143, 'P3', 'P 3', 'P 3', 'P 3', '', '3',
                       'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(431, 144, 'P3_1', 'P 3_1', 'P 3_1', 'P 31', '', '3',
                       'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(432, 145, 'P3_2', 'P 3_2', 'P 3_2', 'P 32', '', '3',
                       'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(433, 146, 'R3', 'R 3', 'R 3', 'R 3', 'H', '3',
                       'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(434, 146, 'R3', 'R 3', 'R 3', 'P 3*', 'R', '3',
                       'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(435, 147, 'P-3', 'P -3', 'P -3', '-P 3', '', '-3',
                       'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(436, 148, 'R-3', 'R -3', 'R -3', '-R 3', 'H', '-3',
                       'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(437, 148, 'R-3', 'R -3', 'R -3', '-P 3*', 'R', '-3',
                       'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(438, 149, 'P312', 'P 3 1 2', 'P 3 1 2', 'P 3 2', '',
                       '32', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(439, 150, 'P321', 'P 3 2 1', 'P 3 2 1', 'P 3 2"', '',
                       '32', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(440, 151, 'P3_112', 'P 3_1 1 2', 'P 3_1 1 2',
                       'P 31 2c (0 0 1)', '', '32', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(441, 152, 'P3_121', 'P 3_1 2 1', 'P 3_1 2 1', 'P 31 2"',
                       '', '32', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(442, 153, 'P3_212', 'P 3_2 1 2', 'P 3_2 1 2',
                       'P 32 2c (0 0 -1)', '', '32', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(443, 154, 'P3_221', 'P 3_2 2 1', 'P 3_2 2 1', 'P 32 2"',
                       '', '32', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(444, 155, 'R32', 'R 3 2', 'R 3 2', 'R 3 2"', 'H', '32',
                       'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(445, 155, 'R32', 'R 3 2', 'R 3 2', 'P 3* 2', 'R', '32',
                       'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(446, 156, 'P3m1', 'P 3 m 1', 'P 3 m 1', 'P 3 -2"', '',
                       '3m', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(447, 157, 'P31m', 'P 3 1 m', 'P 3 1 m', 'P 3 -2', '',
                       '3m', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(448, 158, 'P3c1', 'P 3 c 1', 'P 3 c 1', 'P 3 -2"c', '',
                       '3m', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(449, 159, 'P31c', 'P 3 1 c', 'P 3 1 c', 'P 3 -2c', '',
                       '3m', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(450, 160, 'R3m', 'R 3 m', 'R 3 m', 'R 3 -2"', 'H', '3m',
                       'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(451, 160, 'R3m', 'R 3 m', 'R 3 m', 'P 3* -2', 'R', '3m',
                       'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(452, 161, 'R3c', 'R 3 c', 'R 3 c', 'R 3 -2"c', 'H',
                       '3m', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(453, 161, 'R3c', 'R 3 c', 'R 3 c', 'P 3* -2n', 'R',
                       '3m', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(454, 162, 'P-31m', 'P -3 1 m', 'P -3 1 2/m', '-P 3 2',
                       '', '-3m', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(455, 163, 'P-31c', 'P -3 1 c', 'P -3 1 2/c', '-P 3 2c',
                       '', '-3m', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(456, 164, 'P-3m1', 'P -3 m 1', 'P -3 2/m 1', '-P 3 2"',
                       '', '-3m', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(457, 165, 'P-3c1', 'P -3 c 1', 'P -3 2/c 1', '-P 3 2"c',
                       '', '-3m', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(458, 166, 'R-3m', 'R -3 m', 'R -3 2/m', '-R 3 2"', 'H',
                       '-3m', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(459, 166, 'R-3m', 'R -3 m', 'R -3 2/m', '-P 3* 2', 'R',
                       '-3m', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(460, 167, 'R-3c', 'R -3 c', 'R -3 2/c', '-R 3 2"c', 'H',
                       '-3m', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(461, 167, 'R-3c', 'R -3 c', 'R -3 2/c', '-P 3* 2n', 'R',
                       '-3m', 'trigonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(462, 168, 'P6', 'P 6', 'P 6', 'P 6', '', '6',
                       'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(463, 169, 'P6_1', 'P 6_1', 'P 6_1', 'P 61', '', '6',
                       'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(464, 170, 'P6_5', 'P 6_5', 'P 6_5', 'P 65', '', '6',
                       'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(465, 171, 'P6_2', 'P 6_2', 'P 6_2', 'P 62', '', '6',
                       'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(466, 172, 'P6_4', 'P 6_4', 'P 6_4', 'P 64', '', '6',
                       'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(467, 173, 'P6_3', 'P 6_3', 'P 6_3', 'P 6c', '', '6',
                       'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(468, 174, 'P-6', 'P -6', 'P -6', 'P -6', '', '-6',
                       'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(469, 175, 'P6/m', 'P 6/m', 'P 6/m', '-P 6', '', '6/m',
                       'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(470, 176, 'P6_3/m', 'P 6_3/m', 'P 6_3/m', '-P 6c', '',
                       '6/m', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(471, 177, 'P622', 'P 6 2 2', 'P 6 2 2', 'P 6 2', '',
                       '622', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(472, 178, 'P6_122', 'P 6_1 2 2', 'P 6_1 2 2',
                       'P 61 2 (0 0 -1)', '', '622', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(473, 179, 'P6_522', 'P 6_5 2 2', 'P 6_5 2 2',
                       'P 65 2 (0 0 1)', '', '622', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(474, 180, 'P6_222', 'P 6_2 2 2', 'P 6_2 2 2',
                       'P 62 2c (0 0 1)', '', '622', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(475, 181, 'P6_422', 'P 6_4 2 2', 'P 6_4 2 2',
                       'P 64 2c (0 0 -1)', '', '622', 'hexagonal',
                       'primitive'))
        self.spacegroup.append(
            SpaceGroup(476, 182, 'P6_322', 'P 6_3 2 2', 'P 6_3 2 2', 'P 6c 2c',
                       '', '622', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(477, 183, 'P6mm', 'P 6 m m', 'P 6 m m', 'P 6 -2', '',
                       '6mm', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(478, 184, 'P6cc', 'P 6 c c', 'P 6 c c', 'P 6 -2c', '',
                       '6mm', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(479, 185, 'P6_3cm', 'P 6_3 c m', 'P 6_3 c m', 'P 6c -2',
                       '', '6mm', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(480, 186, 'P6_3mc', 'P 6_3 m c', 'P 6_3 m c',
                       'P 6c -2c', '', '6mm', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(481, 187, 'P-6m2', 'P -6 m 2', 'P -6 m 2', 'P -6 2', '',
                       '-6m2', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(482, 188, 'P-6c2', 'P -6 c 2', 'P -6 c 2', 'P -6c 2',
                       '', '-6m2', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(483, 189, 'P-62m', 'P -6 2 m', 'P -6 2 m', 'P -6 -2',
                       '', '-6m2', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(484, 190, 'P-62c', 'P -6 2 c', 'P -6 2 c', 'P -6c -2c',
                       '', '-6m2', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(485, 191, 'P6/mmm', 'P 6/m m m', 'P 6/m 2/m 2/m',
                       '-P 6 2', '', '6/mmm', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(486, 192, 'P6/mcc', 'P 6/m c c', 'P 6/m 2/c 2/c',
                       '-P 6 2c', '', '6/mmm', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(487, 193, 'P6_3/mcm', 'P 6_3/m c m', 'P 6_3/m 2/c 2/m',
                       '-P 6c 2', '', '6/mmm', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(488, 194, 'P6_3/mmc', 'P 6_3/m m c', 'P 6_3/m 2/m 2/c',
                       '-P 6c 2c', '', '6/mmm', 'hexagonal', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(489, 195, 'P23', 'P 2 3', 'P 2 3', 'P 2 2 3', '', '23',
                       'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(490, 196, 'F23', 'F 2 3', 'F 2 3', 'F 2 2 3', '', '23',
                       'cubic', 'face'))
        self.spacegroup.append(
            SpaceGroup(491, 197, 'I23', 'I 2 3', 'I 2 3', 'I 2 2 3', '', '23',
                       'cubic', 'body'))
        self.spacegroup.append(
            SpaceGroup(492, 198, 'P2_13', 'P 2_1 3', 'P 2_1 3', 'P 2ac 2ab 3',
                       '', '23', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(493, 199, 'I2_13', 'I 2_1 3', 'I 2_1 3', 'I 2b 2c 3',
                       '', '23', 'cubic', 'body'))
        self.spacegroup.append(
            SpaceGroup(494, 200, 'Pm3', 'P m 3', 'P 2/m -3', '-P 2 2 3', '',
                       'm-3', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(495, 201, 'Pn3', 'P n 3', 'P 2/n -3', 'P 2 2 3 -1n',
                       '1', 'm-3', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(496, 201, 'Pn3', 'P n 3', 'P 2/n -3', '-P 2ab 2bc 3',
                       '2', 'm-3', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(497, 202, 'Fm3', 'F m 3', 'F 2/m -3', '-F 2 2 3', '',
                       'm-3', 'cubic', 'face'))
        self.spacegroup.append(
            SpaceGroup(498, 203, 'Fd3', 'F d 3', 'F 2/d -3', 'F 2 2 3 -1d',
                       '1', 'm-3', 'cubic', 'face'))
        self.spacegroup.append(
            SpaceGroup(499, 203, 'Fd3', 'F d 3', 'F 2/d -3', '-F 2uv 2vw 3',
                       '2', 'm-3', 'cubic', 'face'))
        self.spacegroup.append(
            SpaceGroup(500, 204, 'Im3', 'I m 3', 'I 2/m -3', '-I 2 2 3', '',
                       'm-3', 'cubic', 'body'))
        self.spacegroup.append(
            SpaceGroup(501, 205, 'Pa3', 'P a 3', 'P 2_1/a -3', '-P 2ac 2ab 3',
                       '', 'm-3', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(502, 206, 'Ia3', 'I a 3', 'I 2_1/a -3', '-I 2b 2c 3',
                       '', 'm-3', 'cubic', 'body'))
        self.spacegroup.append(
            SpaceGroup(503, 207, 'P432', 'P 4 3 2', 'P 4 3 2', 'P 4 2 3', '',
                       '432', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(504, 208, 'P4_232', 'P 4_2 3 2', 'P 4_2 3 2',
                       'P 4n 2 3', '', '432', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(505, 209, 'F432', 'F 4 3 2', 'F 4 3 2', 'F 4 2 3', '',
                       '432', 'cubic', 'face'))
        self.spacegroup.append(
            SpaceGroup(506, 210, 'F4_132', 'F 4_1 3 2', 'F 4_1 3 2',
                       'F 4d 2 3', '', '432', 'cubic', 'face'))
        self.spacegroup.append(
            SpaceGroup(507, 211, 'I432', 'I 4 3 2', 'I 4 3 2', 'I 4 2 3', '',
                       '432', 'cubic', 'body'))
        self.spacegroup.append(
            SpaceGroup(508, 212, 'P4_332', 'P 4_3 3 2', 'P 4_3 3 2',
                       'P 4acd 2ab 3', '', '432', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(509, 213, 'P4_132', 'P 4_1 3 2', 'P 4_1 3 2',
                       'P 4bd 2ab 3', '', '432', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(510, 214, 'I4_132', 'I 4_1 3 2', 'I 4_1 3 2',
                       'I 4bd 2c 3', '', '432', 'cubic', 'body'))
        self.spacegroup.append(
            SpaceGroup(511, 215, 'P-43m', 'P -4 3 m', 'P -4 3 m', 'P -4 2 3',
                       '', '-43m', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(512, 216, 'F-43m', 'F -4 3 m', 'F -4 3 m', 'F -4 2 3',
                       '', '-43m', 'cubic', 'face'))
        self.spacegroup.append(
            SpaceGroup(513, 217, 'I-43m', 'I -4 3 m', 'I -4 3 m', 'I -4 2 3',
                       '', '-43m', 'cubic', 'body'))
        self.spacegroup.append(
            SpaceGroup(514, 218, 'P-43n', 'P -4 3 n', 'P -4 3 n', 'P -4n 2 3',
                       '', '-43m', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(515, 219, 'F-43c', 'F -4 3 c', 'F -4 3 c', 'F -4c 2 3',
                       '', '-43m', 'cubic', 'face'))
        self.spacegroup.append(
            SpaceGroup(516, 220, 'I-43d', 'I -4 3 d', 'I -4 3 d',
                       'I -4bd 2c 3', '', '-43m', 'cubic', 'body'))
        self.spacegroup.append(
            SpaceGroup(517, 221, 'Pm-3m', 'P m -3 m', 'P 4/m -3 2/m',
                       '-P 4 2 3', '', 'm-3m', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(518, 222, 'Pn-3n', 'P n -3 n', 'P 4/n -3 2/n',
                       'P 4 2 3 -1n', '1', 'm-3m', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(519, 222, 'Pn-3n', 'P n -3 n', 'P 4/n -3 2/n',
                       '-P 4a 2bc 3', '2', 'm-3m', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(520, 223, 'Pm-3n', 'P m -3 n', 'P 4_2/m -3 2/n',
                       '-P 4n 2 3', '', 'm-3m', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(521, 224, 'Pn-3m', 'P n -3 m', 'P 4_2/n -3 2/m',
                       'P 4n 2 3 -1n', '1', 'm-3m', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(522, 224, 'Pn-3m', 'P n -3 m', 'P 4_2/n -3 2/m',
                       '-P 4bc 2bc 3', '2', 'm-3m', 'cubic', 'primitive'))
        self.spacegroup.append(
            SpaceGroup(523, 225, 'Fm-3m', 'F m -3 m', 'F 4/m -3 2/m',
                       '-F 4 2 3', '', 'm-3m', 'cubic', 'face'))
        self.spacegroup.append(
            SpaceGroup(524, 226, 'Fm-3c', 'F m -3 c', 'F 4/m -3 2/c',
                       '-F 4c 2 3', '', 'm-3m', 'cubic', 'face'))
        self.spacegroup.append(
            SpaceGroup(525, 227, 'Fd-3m', 'F d -3 m', 'F 4_1/d -3 2/m',
                       'F 4d 2 3 -1d', '1', 'm-3m', 'cubic', 'face'))
        self.spacegroup.append(
            SpaceGroup(526, 227, 'Fd-3m', 'F d -3 m', 'F 4_1/d -3 2/m',
                       '-F 4vw 2vw 3', '2', 'm-3m', 'cubic', 'face'))
        self.spacegroup.append(
            SpaceGroup(527, 228, 'Fd-3c', 'F d -3 c', 'F 4_1/d -3 2/c',
                       'F 4d 2 3 -1cd', '1', 'm-3m', 'cubic', 'face'))
        self.spacegroup.append(
            SpaceGroup(528, 228, 'Fd-3c', 'F d -3 c', 'F 4_1/d -3 2/c',
                       '-F 4cvw 2vw 3', '2', 'm-3m', 'cubic', 'face'))
        self.spacegroup.append(
            SpaceGroup(529, 229, 'Im-3m', 'I m -3 m', 'I 4/m -3 2/m',
                       '-I 4 2 3', '', 'm-3m', 'cubic', 'body'))
        self.spacegroup.append(
            SpaceGroup(530, 230, 'Ia-3d', 'I a -3 d', 'I 4_1/a -3 2/d',
                       '-I 4bd 2c 3', '', 'm-3m', 'cubic', 'body'))

        # temporary
        for g in self.spacegroup:
            # trigonal rhombohedral
            # hexagonal trigonal
            if(g.number == 143 or g.number == 144 or g.number == 145 \
                    or g.number == 147 or g.number == 149 or g.number == 150 \
                    or g.number == 151 or g.number == 152 or g.number == 153 \
                    or g.number == 154 or g.number == 156 or g.number == 157 \
                    or g.number == 158 or g.number == 159 or g.number == 162 \
                    or g.number == 163 or g.number == 164 or g.number == 165):
                g.shape = 'hexagonal'
            if (g.center == 'primitive'):
                g.center = 'simple'

    def findByName(self, name):
        for g in self.spacegroup:
            name0 = g.name
            name0 = name0.replace(' ', '')
            name0 = name0.replace('_', '')
            if (name0 == name):
                return g
        for g in self.spacegroup:
            name0 = g.name
            name0 = name0.replace(' ', '')
            name0 = name0.replace('_', '')
            if (g.nameHM == name):
                return g
        for g in self.spacegroup:
            name0 = g.name
            name0 = name0.replace(' ', '')
            name0 = name0.replace('_', '')
            if (g.nameFull == name):
                return g
        return SpaceGroup()

    def findByHallNumber(self, number):
        for g in self.spacegroup:
            if (g.hallnumber == number):
                return g
        return SpaceGroup()

    def findByIntTableNumber(self, number):
        for g in self.spacegroup:
            if (g.number == number):
                return g
        return SpaceGroup()


if __name__ == '__main__':
    spgtable = SpaceGroupTable()
    print(
        'hallnumber',
        'number',
    )
    print(
        'name',
        'nameHM',
        'nameHM',
        'nameFull',
        'nameHall',
        'axis_code',
        'pointgroup',
    )
    print('shape', 'center')
    for sg in spgtable.spacegroup:
        print(sg.hallnumber, )
        print(sg.number, )
        print(sg.name, )
        print(sg.nameHM, )
        print(sg.nameFull, )
        print(sg.nameHall, )
        print(sg.axis_code, )
        print(sg.pointgroup, )
        print(sg.shape, )
        print(sg.center)
