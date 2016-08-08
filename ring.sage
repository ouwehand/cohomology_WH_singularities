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


#This file contains the definition of the polynomial ring, which is the only global definition
#
#The number of variables and the variable names may be changed. The monomial order can be changed to any other global monomial order
#
#The name 'R' should not be modified
#
#The base field must remain \Q because SAGE has very limited support for polynomials over p-adic fields

#Note: the monomial order 'degrevlex' is probably a better choice than 'lex'. Also see the file order.sage. 
#
#But for now we will stick to 'lex' because that is what appears in the algorithm listings in paragraphs 4.2.3 and 4.2.4

R.<x,y,z,u> = PolynomialRing(QQ, order='lex')

vars = R.gens()

