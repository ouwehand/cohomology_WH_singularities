
About
-------

This SAGE program approximates the Frobenius action on the `G(\underline{w})`-invariant part of the rigid cohomology space `H^{n-1}(\mathbb{P}^{n-1} \setminus \tilde{S}_{\infty})`.

It is based on a modified version of the algorithm by Abbott, Kedlaya and Roe, see _Bounding Picard numbers of surfaces using p-adic cohomology_.

Coded by David Ouwehand. This program is released under the terms of the GNU General Public License.

The code is meant to be easy to understand and easy to modify. It is not meant to be as efficient as possible. The main purpose is to carry out simple experiments that cannot (yet) be computed by other implementations; see the file `tests.sage`.

To get started quickly, start up SAGE in the same directory as the source files and type: 

`load("tests.sage")`

You can also try out this code on the SAGE cloud: <https://cloud.sagemath.com/>.


References
----------

The algorithm is explained in detail in the author's Ph.D. thesis: _Local rigid cohomology of weighted homogeneous hypersurface singularities_.

The most relevant parts are section 3.1 (basic definitions, also outlined below) and section 4.2 (statement and proof of the algorithm).

On a theoretical level, the code given here is (more or less) equivalent to the _Frobenius Project_ by Aise Johan de Jong:

<http://math.columbia.edu/~dejong/algebraic_geometry/Frobenius/>


Definitions / Inputs
--------------------

A weighted homogeneous polynomial G with integer coefficients.

The variables of G are written: `vars = [x_1, ..., x_n]`.

The weights of G are written: `weights = [w_1, ..., w_n]`.

The weighted degree of G is written d.

Fix a prime number p. Also consider the finite field `F_p` and the p-adic field `Q_p`.

The projective hypersurface `\tilde{S}_{\infty}` over `F_p` is defined by the equation `G(x_1^{w_1}, ..., x_n^{w_n})` modulo p.

The group `G(\underline{w})` is defined as the product `(Z / w_1 Z) x ... x (Z / w_n Z)`. The rigid cohomology space `H^{n-1}(\mathbb{P}^{n-1} \setminus \tilde{S}_{\infty})` comes with an action of this group.


Assumptions
-----------

The weighted degree d is not divisible by p.

None of the weights `w_i` is divisible by p.

The projective hypersurface `\tilde{S}_{\infty}` over `F_p` is smooth.

The projective hypersurface over `Q_p`, defined by the equation `G(x_1^{w_1}, ..., x_n^{w_n})`, is smooth.


Calculation / Output
--------------------

To approximate the Frobenius matrix on the `G(\underline{w})`-invariant cohomology, use the function

`frobenius_matrix(G, weights, p, N, digits)`.

The parameters G, weights, p are as above.

N is the index at which to truncate the Frobenius action on overconvergent differentials.

`digits` is the number of p-adic digits to use in the final answer.


Overview of source files
------------------------

frobenius.sage: Contains the function `frobenius_matrix` described above.

basis.sage: Computes a basis for the `G(\underline{w})`-invariant cohomology.

ring.sage: Only contains the definition of the polynomial ring. Here you can change the number of variables and the monomial order.

tests.sage: Contains a few test cases / examples.

dimension.sage: Computes the dimension of the `G(\underline{w})`-invariant cohomology. This is completely independent from all the other source files.

order.sage: This source file is only meant to demonstrate the remarks in paragraph 4.2.5 of the author's thesis.


Comparison with the Frobenius Project
-------------------------------------

This program was originally created to verify that the algorithms in section 4.2 of the author's thesis are really equivalent to the algorithm behind the Frobenius Project (see in particular paragraph 4.2.6 for details). That is why the implementation was kept as basic as possible. In the end, this implementation has the advantage that it works under a *different* set of assumptions/restrictions.

Here are the assumptions of the implementation by de Jong:

* The number of variables is equal to 3 or 4.

* The weighted degree d is not smaller than the sum of the weights.

* The monomial order is either the lexicographic one or the graded reverse lexicographic one.

For examples that satisfy the assumptions of *both* implementations, the version of de Jong will compute the answer more quickly.

