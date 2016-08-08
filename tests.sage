#    Copyright (C) 2016  David Ouwehand
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>. 


#A few examples / test cases

load("ring.sage")
load("dimension.sage")
load("basis.sage")
load("frobenius.sage")

##################################################################################

#Test 1: Griffiths-Dwork reduction in the homogeneous case
#
#This example is simple enough to verify by hand

print("test1")
print("")

G = x^3 + y^3 + z^3 + u^3

v = basis_below(G, [1, 1, 1, 1])

print( reduce_pole_order( x*y*G + x^4*y, 3, v[0], v[1], v[2] ) )

#Answer: [0, 0, 0, 0, 0, 4/3]


##################################################################################

#Test 2: Local cohomology of an ordinary double point
#
#See proposition 4.3.6

print("")
print("test2")
print("")

G = x^2 + y^2 + z^2 + u^2

print( frobenius_matrix( G, [1, 1, 1, 1], 3, 5, 5 ) )

#Answer: 9 \in \Q_3


##################################################################################

#Test 3: Local cohomology of a singularity of type A_3
#
#See proposition 4.3.10

print("")
print("test3")
print("")

G = x^4 + y^2 + z^2 + u^2

print( frobenius_matrix( G, [1, 2, 2, 2], 3, 5, 5 ) )

#Answer: 9 \in \Q_3


##################################################################################

#Test 4: An example in two variables
#
#See proposition 4.3.12 and remark 4.3.14

print("")
print("test4")
print("")

R.<x,y> = PolynomialRing(QQ, order='lex')
vars = R.gens()

G = x^3 + y^3

print( frobenius_matrix( G, [1, 1], 5, 10, 5 ).charpoly('t') )

#Answer: t^2 - 5^2

print( frobenius_matrix( G, [1, 1], 7, 10, 5 ).charpoly('t') )

#Answer: (t - 7)^2


##################################################################################

#The next two tests should be compared with the default example in the Frobenius Project by Aise Johan de Jong
#
#http://math.columbia.edu/~dejong/algebraic_geometry/Frobenius/
#
#Also see paragraph 4.2.6

R.<x,y,z,u> = PolynomialRing(QQ, order='lex')
vars = R.gens()

print("")
print("test5")
print("")

print( compute_dimension( 64, [8, 7, 15, 19] ) )

#Answer: 7

print("")
print("test6")
print("")

#Warning: \tilde{G} is singular!

G = -x^8 + x^3*y^3*u + x*y^8 - x*y*z^2*u - y^7*z + y*u^3 + z^3*u

v = basis_below(G, [8, 7, 15, 19])

print(v[2])

#Answer: [z, x*y, y*z*u^3, z^4*u, y^7*z^2, y*z^4*u^4, z^7*u^2]


##################################################################################

#Test 7: A weighted Dwork hypersurface

print("")
print("test7")
print("")

R.<x,y,z> = PolynomialRing(QQ, order='lex')
vars = R.gens()

G = x^5 + y^10 + z^2 + x*y^3*z

print( frobenius_matrix( G, [2, 1, 5], 3, 10, 5 ).charpoly('t') )

# The approximate characteristic polynomial seems to converge to:
#
# t^4 - 3 * t^3 - 3^4 * t + 3^6

# The Frobenius matrix above can also be compared with the matrix that is computed by the 'curves' variation of the Frobenius Project
#
# But be aware that the matrix computed by de Jong's code has an extra factor 3^{-1} in front of it
 
