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


load("ring.sage")
load("basis.sage")

### Auxiliary functions ###

#Add the array 'b' to the array 'a'
def add_to_array(b, a):
	
	for (i, elt) in enumerate(b):
		a[i] += elt
	return

#The q-power Frobenius on polynomials
def frobenius(A, q):
	
	return sum( c * m^q for (c, m) in zip( A.coefficients(), A.monomials() ) )


#Computes the cohomology class of the form ( h * \tilde{A} * Omega ) / ( \tilde{G}^j ) as a linear combination of the basis, where
#
# h = prod_i x_i^{w_i - 1} 
#
#The argument 'I' is a list containing the partial derivatives of G, *not* of \tilde{G}
#
#The argument 'Q' is a Groebner basis for 'I'.
#
#The arguments 'I', 'Q' and 'basis' should be computed by 'basis_below'. This ensures that the remainder 'A rem Q' from the multivariate division algorithm is a linear combination of the elements of 'basis'.
#
#The function sends back a list 'coeffs' where coeffs[i] is the coefficient for basis[i]

def reduce_pole_order(A, j, I, Q, basis):

	coeffs = [0 for b in basis]
	
	#Compute the remainder
	r = A.reduce(Q)

	for (m, c) in zip( r.monomials(), r.coefficients() ):
		coeffs[ basis.index( m ) ] = c
	
	if A != r:

		#Write A - r as a combination of the partial derivatives
		quotients = (A-r).lift(I)
	
		#Compute partial derivatives of the quotients
		B = sum( derivative(quot, var) for (quot, var) in zip( quotients, vars ) )

		#Apply the rule for pole order reduction
		add_to_array( reduce_pole_order( (1 / (j-1) ) * B, j-1, I, Q, basis ), coeffs )
	
	return coeffs


#An iterative implementation of the pole order reduction function
#
#This version should be slightly faster. It doesn't follow the structure from the text as closely as the previous version.
#
#The algorithm has been broken up into two parts. 
#In this way the code for pole orders j > n-1 (where the bulk of the computation happens) is as fast as possible.

def reduce_pole_order_iter(A, j, I, Q, basis):

	coeffs = [0 for b in basis]
	n = len(vars)

	#If j > n-1 then we only need to reduce the pole order
	
	for j in range(j, n-1, -1):

		quotients = A.lift(I)

		A = (1 / (j-1)) * sum( derivative(quot, var) for (quot, var) in zip( quotients, vars ) )

	#Only forms of pole order <= n-1 contribute to the coefficient vector

	for j in range( min(n-1, j), 0, -1):
	
		#Compute the remainder
		r = A.reduce(Q)

		for (m, c) in zip( r.monomials(), r.coefficients() ):
			coeffs[ basis.index( m ) ] += c
		
		#This test is stronger than 'if j==1': the basis may run out before the pole order reaches the value j=1
		if A == r:
			break

		#Write A - r as a combination of the partial derivatives
		quotients = (A-r).lift(I)
		
		#Reduce the pole order
		A = (1 / (j-1)) * sum( derivative(quot, var) for (quot, var) in zip( quotients, vars ) )
	
	return coeffs


#Compute the Frobenius matrix.
#
#Only works for q = p for a prime p
#
#N is the index at which to truncate the Frobenius action on overconvergent differentials.
#
#'digits' is the number of p-adic digits to use in the final answer.

#This function admits some obvious optimizations
#However, these are not worth it because the bulk of the computation happens in the pole order reduction step

def frobenius_matrix(G, weights, p, N, digits):

	#Compute the basis
	
	v = basis_below(G, weights)
	I = v[0]
	Q = v[1]
	basis = v[2]

	#Some useful quantities

	s = sum(weights)
	d = weighted_degree(G.monomials()[0], weights)

	#Define the matrix (with entries in a p-adic field) that will hold the answer
	
	K = Qp( p, digits)
	M = matrix(K, len(basis), 0)

	#'column' will contain the next column of 'M' before we copy it over
	column = [0 for b in basis]

	#This array will contain parts of the numerators of the terms in the truncated sum that approximates the Frobenius

	numerators = [0 for i in range(0, N) ]

	#Start by defining 'numerators'. This is independent of the basis elements
	
	prefix = p^(len(vars) - 1) *  prod( var^(p-1) for var in vars )
	numerators[0] = prefix

	pDelta = G^p - frobenius(G, p)

	for i in range(1, N):

		numerators[i] = numerators[i-1] * pDelta

	#Start to loop over the basis elements. 

	for m in basis:
		
		#Determine the pole order corresponding to 'm'
		t = (weighted_degree(m, weights) + s ) / d

		#Now compute the total numerator of the truncated sum that approximates the Frobenius
		total_numerator = 0

		for (i, num) in enumerate( numerators ):
			
			total_numerator = G^p * total_numerator + num * binomial(t+i-1, t-1)

		total_numerator *= frobenius(m, p)

		#Reduce the pole order and add the result as a column to 'M'
		
		column = reduce_pole_order_iter( total_numerator, p * (N - 1 + t), I, Q, basis ) 

		M = M.augment( vector(K, column) )

	return	M

