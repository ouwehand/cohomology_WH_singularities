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


#Compute the dimension of local cohomology of a singularity that is given by a polynomial of degree 'd' with respect to weights 'weights'.
#
#This implementation assumes that the length of the list 'weights' is at least 3. When 'weights' has length 2 the dimension of local cohomology is one larger than the computed value.
#
#However, the value computed by the function below is always equal to the length of the basis that is computed by basis_below.

def compute_dimension(d, weights):
	
	n = len(weights)
	s = sum(weights)

	#It may happen that the first entry of the 'max' is zero. Then we set 'precision' to 1 in order to avoid a crash
	precision = max((n-1) * d - s + 1, 1)
	R.<t> = PowerSeriesRing(ZZ, default_prec = precision) 
	
	Psi = prod( (t^(d - w) - 1) / (t^w - 1) for w in weights )

	powers_of_t = ( i * d - s for i in range(1, n) if i * d - s >= 0 )

	return sum( Psi[p] for p in powers_of_t )


