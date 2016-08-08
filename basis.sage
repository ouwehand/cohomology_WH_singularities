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

### Auxiliary functions ###


#Apply the tilde operator to a monomial.
def tilde(monomial, weights):
	
	return prod( var^(exponent * weight) for (var, exponent, weight) in zip(vars, monomial.exponents()[0], weights) )

#Compute the weighted degree of a monomial
def weighted_degree(monomial, weights):

	return sum( exponent * weight for (exponent, weight) in zip(monomial.exponents()[0], weights) )

#Check if a monomial is divisible by some element in a list of monomials
def divisible_by_element(monomial, lst_monomials):
	
	for m in lst_monomials:
		
		#Temporarily set this flag to 'True'
		m_is_divisor = True

		for (exp1, exp2) in zip( monomial.exponents()[0], m.exponents()[0] ):

			if exp1 < exp2:
				
				#The current candidate 'm' is not a divisor of 'monomial'
				m_is_divisor = False
				break

		if m_is_divisor == True:
			return True
	
	#If the main loop ends then we have found no divisor
	return False


### Implementation of the basis algorithm ###

#Returns a list 'v' where
#
# v[0] are the partial derivatives of 'G'
#
# v[1] is the reduced Groebner basis for the Jacobian ideal of G
#
# v[2] is the basis "below".

def basis_below(G, weights):
	
	d = weighted_degree(G.monomials()[0], weights)
	s = sum(weights)
	n = len(vars)
	derivatives = [derivative(G, var) for var in vars]
	b = prod(vars)
	basis = []

	#Compute the Groebner basis of the Jacobian ideal of G
	Q = R.ideal(derivatives).groebner_basis()
	leading_terms = [q.lt() for q in Q]

	#Compute the degrees that correspond to the various pole orders
	
	degrees = []
	for i in range(0, n):
		D = i*d - s

		if D > 0:
			degrees.append(D)
		if D == 0:
			#For D = 0 we add the element '1' to the basis 
			basis.append( R(1) )

	for D in degrees:
		
		#Compute the set of monomials of weighted degree 'D' that are not  a multiple of an element in 'leading_terms'
		#
		#This set is added to the basis

		exponents = WeightedIntegerVectors(D, weights)
		monomials = (tilde(b, e) for e in exponents)

		basis += [ m for m in monomials if not divisible_by_element(m, leading_terms) ]

	return [ derivatives, Q, basis ]

		

			

