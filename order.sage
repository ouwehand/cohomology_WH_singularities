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


#This source file demonstrates what happens when the monomial order is changed
#
#It is only meant to demonstrate the remarks in paragraph 4.2.5 of the author's thesis

load("basis.sage")
load("frobenius.sage")

### Auxiliary functions ###


#Compute the tilde of a polynomial, rather than a monomial
def tilde_poly(A, weights):
	
	return sum( c * tilde(m, weights) for (c,m) in A )


#Check if the condition from definition 4.2.1 holds for a given monomial
def check_degrees_condition(m, weights):

	flag = True

	for (d,w) in zip( m.degrees(), weights):
		
		if mod(d, w) != w-1:
			flag = false
			break
	return flag


### Definitions of ideals and Groebner bases ###

#The ring 'R' and the equation 'G' may be modified
#
#This will not alter the functionality of the program, except for two places where manual changes are required


#Definition of the polynomial ring

R.<x,y,z> = PolynomialRing(QQ, order='deglex')
vars = R.gens()

#Definition of the equation and weights

G = x^5 + y^10 + z^2 + x*y^3*z
weights = [2, 1, 5]
prime = 3

#The Groebner bases of the Jacobian ideals I and J

v = basis_below(G, weights)

Gt = tilde_poly(G, weights)
vv = basis_below(Gt , [1 for var in vars] )

I = R.ideal( v[0] )
J = R.ideal( vv[0] )

QI = v[1]
QJ = vv[1]

#Definition of the ideal \tilde{I}

It = R.ideal( [ tilde_poly(der, weights) for der in I.gens() ] )
QIt = It.groebner_basis()

#Compute the ideals of leading terms

LI = R.ideal( [ q.lt() for q in QI ] )
LIt = R.ideal( [ q.lt() for q in QIt ] )
LJ = R.ideal( [ q.lt() for q in QJ ] )

#Compute the monomials \tilde{u}_j

Ut = [ tilde(u, weights) for u in LI.gens() ]

##################################################################################

#Test 1 shows that proposition 4.2.9 breaks down

print("\n\n test 1: Groebner basis of It \n")

print( [ tilde_poly(q, weights) for q in QI ] )
print("")
print( QIt )

##################################################################################

# Test 2: the proof of corollary 4.2.10 breaks down
#
# Ut \not \subset LIt
#
# Which implies:
#
# { \tilde{m} | m is a leading monomial of I} \not \subset LIt

print("\n\n test 2: leading term ideals of I and It \n")

print("Ut: \n")
print(Ut)
print("")

print("LIt: \n")
print(LIt)
print("")

print( [ (ut in LIt) for ut in Ut ] )

##################################################################################

# Test 3: the statement of proposition 4.2.8 is no longer true
#
# If u belongs to LI, then \varphi(u) does not necessarily belong to LJ
#
# This shows that the image of \mathcal{N}_D under \varphi is not equal to \mathcal{M}'_D

print("\n\n test 3: The image of varphi \n")

b = prod( vars )
h = tilde( b, [w-1 for w in weights] )

#Note: we only test the generators of LI
#
#Of course the ideal LI contains many more monomials
#
#The value 'False' appears in the list below
#
#This means that \mathcal{M}'_D is generally not a subset of \varphi(\mathcal{N}_D)

print("Generators of LI: \n")
print(LI.gens())
print("")

print( [ (h * tilde(u, weights) in LJ) for u in LI.gens() ] )
print("")

#The following test shows that \varphi(\mathcal{N}_D) is generally not a subset of \mathcal{M}'_D
#
#NOTE: The lines below only work for the default example

print("Generators of LJ: \n")
print(LJ.gens())
print("")

mon1 = x^3*y^2*z^9
mon2 = x*y^2*z

print( mon1 in LJ )
print( mon1 == h * tilde( mon2 , weights) )
print( mon2 in LI )

##################################################################################

# Test 4: Comparison of the canonical basis (defined by the sets \mathcal{M}'_D) and the basis computed by the algorithm

print("\n\n test 4: Computed basis vs canonical basis \n")

#Compute the canonical basis using the naive method (pruning the larger basis upstairs)

basis_homogeneous = vv[2]

basis_canonical = [ m for m in basis_homogeneous if check_degrees_condition(m, weights) ]

#The basis computed by the algorithm

basis_algorithm = [ h * tilde(m, weights) for m in v[2] ]

print("Canonical basis: \n")
print(basis_canonical)
print("")

print("Computed basis: \n")
print(basis_algorithm)
print("")

#see if the two bases are the same

flag = True

copy_can = [m for m in basis_canonical]
copy_can.sort()

copy_algo = [m for m in basis_algorithm]
copy_algo.sort()

if len(copy_can) != len(copy_algo):
	flag = False
else:
	for (m,n) in zip(copy_can, copy_algo):

		if m != n:
			flag = False
			break

print(flag)

##################################################################################

#Test 5: See if basis_algorithm is a basis for the G(w)-invariant cohomology

print("\n\n test 5: Check that basis_algorithm is a basis for the G(w)-invariant cohomology \n")

#See if the two sets have the same number of elements

print( len(basis_algorithm) == len(basis_canonical) )
print("")

#Write each element of basis_algorithm as a linear combination of basis_canonical

M = matrix(QQ, len(basis_canonical), 0)

for m in basis_algorithm:

	ord = ( m.degree() + len(vars) ) / Gt.degree()
	
	coords_homogeneous = reduce_pole_order_iter( m, ord, vv[0], vv[1], vv[2] )

	coords_canonical = [ c for (c, m) in zip(coords_homogeneous, basis_homogeneous) if m in basis_canonical ]

	M = M.augment( vector(QQ, coords_canonical) )

#Check that the resulting matrix is invertible
#
#This means that basis_algorithm is indeed a basis

print( M.determinant() )

##################################################################################

# Test 6: Griffiths-Dwork reduction w.r.t. basis_algorithm

print("\n\n test 6: Griffiths-Dwork reduction w.r.t. basis_algorithm \n")

for ord in range(1, 5):

	#Generate a random polynomial corrpesponding to pole order 'ord'

	d = ( ord * Gt.degree() ) - sum( weights )

	degrees = WeightedIntegerVectors( d, weights )

	poly = sum( randint(-1, 1) * tilde(b, degreevector) for degreevector in degrees )

	#Use modified Griffith-Dwork "below" 
	#
	#This should calculate the coordinates w.r.t. basis_algorithm

	coords_algorithm = reduce_pole_order_iter( poly, ord, v[0], v[1], v[2] )

	#Use classical Griffiths-Dwork "above" to calculate the coordinates w.r.t. basis_canonical

	coords_homogeneous = reduce_pole_order_iter( h * tilde_poly( poly, weights), ord, vv[0], vv[1], vv[2] )

	coords_canonical = [ c for (c, m) in zip(coords_homogeneous, basis_homogeneous) if m in basis_canonical ]

	#Use the base-change matrix from test 5 to verify that 'coords_algorithm' are indeed the right coordinates

	vect1 = vector(QQ, coords_canonical).column()

	vect2 = vector(QQ, coords_algorithm).column()

	print( vect1 == M * vect2 )


##################################################################################

# Test 7: Comparison of the characteristic polynomials of Frobenius for different monomial orders

print("\n\n test 7: The characteristic polynomial of Frobenius \n")

N = 10

P = frobenius_matrix(G, weights, prime, N, 10).characteristic_polynomial('t')

print( P )
print("")

#Compare this approximate polynomial with the one that you get for order='lex'
#
#NOTE: The lines below only work for the default example

R.<x,y,z> = PolynomialRing(QQ, order='lex')
vars = R.gens()
G = x^5 + y^10 + z^2 + x*y^3*z

Plex = frobenius_matrix(G, weights, prime, N, 10).characteristic_polynomial('t')
print( Plex )

#For the default example (prime=3), both monomial orders produce the polynomial
#
# t^4 - prime * t^3 - prime^4 * t + prime^6

#But observe that the calculations with 'deglex' give more precision than 'lex', for the same level of truncation N
#
#The same thing happens if you compare 'degrevlex' with 'lex', but in this case basis_algorithm and basis_canonical are the same!
