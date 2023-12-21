# cliffalg

## Multivector data type for python 3.8

Contains a class for real Clifford algebras with arbitrary basis vectors and another class for multivectors in a specified Clifford algebra.

### Blade naming

The names for the basis blades in an algebra are generated with its constituent vectors sorted alphabetically **not** cyclically permutated. I.e. $xy$, $xz$, $yz$, **not** $xy$, $yz$, $zx$.

### Arithmetic operations

+ Supports addition (+), subtraction (-) and Clifford product (*) with other multivectors of the same algebra or with real scalars i.e. ints or floats.

+ A multivector may be divided (/) by a real scalar but no type may be divided by a multivector.

+ Negation is supported.

### Exponentiation

+ Nonnegative integer powers of multivectors are supported. No types may be exponentiated by a multivector using the **-operator.

+ A method exists for exp of a multivector. It is calculated via the Taylor series for $e^{x}$. If a multivector $\boldsymbol{A}$ is ofthe fom $\boldsymbol{A} = t\boldsymbol{J}$, where $t$ is a positive scalar and $\boldsymbol{J}^{2} = -1$, then it is instead calculated: $\exp(\boldsymbol{A}) = \cos(t) + \boldsymbol{J}\sin(t)$. This is useful for constructing rotors.

## Basic theory

A multivector is any element in a Clifford algebra. A Clifford algebra may be constructed from a vector space $V = \text{span}\{\boldsymbol{e_{1}},\boldsymbol{e_{2}},\dots,\boldsymbol{e_{n}}\}$ by assigning a number to each basis vector $\boldsymbol{e_{i}}$, called its signature. The Clifford product is then defined such that the product of a basis vector with itself is its signature.

The commonly studied Clifford ALgebras only have $1$ or $-1$ as signatures. They are usually written $\mathcal{Cl}_{p,\,n}(\mathbb{F})$, where $p$ is the number of basis vectors with signature $1$, $n$ is the number with signature $-1$ and $\mathbb{F}$ is the field over which the vector space spans.

The product of two different basis vectors is anticommutative ($\boldsymbol{e_{1}e_{2}} = -\boldsymbol{e_{2}e_{1}}$).

The product of $k$ linearly independent vectors is called a $k$-blade and we call $k$ the *grade* of the blade. A $0$-blade is a scalar in $\mathbb{F}$. Multivectors are linear combinations of all possible blades with grades $0\le k \le n$. As an example, a general multivector in the *space algebra* $\mathcal{Cl}_{3,\,0}(\mathbb{R})$ can be written

$$ \boldsymbol{A} = c_{0} + c_{1}\boldsymbol{e_{1}} + c_{2}\boldsymbol{e_{2}} + c_{3}\boldsymbol{e_{3}} + c_{12}\boldsymbol{e_{1}e_{2}} + c_{23}\boldsymbol{e_{2}e_{3}} + c_{31}\boldsymbol{e_{3}e_{1}} + c_{123}\boldsymbol{e_{1}e_{2}e_{3}} $$

where the coefficients are real numbers. The first term is called the scalar part, the following three the vector part, the next three are called the *bivector* part and the last term is the *trivector* part.

The Clifford product of two multivectors is easily calculated using the signatures and the anti-commutation rule. As an example, let $\boldsymbol{A} = \boldsymbol{e_{1}}+\boldsymbol{e_{2}}$, $\boldsymbol{B} = \boldsymbol{e_{1}}+\boldsymbol{e_{2}}+\boldsymbol{e_{3}}$ be elements of $\mathcal{Cl}_{3,\,0}(\mathbb{R})$. Then

$$\begin{aligned}\boldsymbol{AB}
&= (\boldsymbol{e_{1}} + \boldsymbol{e_{2}})(\boldsymbol{e_{1}} + \boldsymbol{e_{2}} + \boldsymbol{e_{3}}) 
= \boldsymbol{e_{1}e_{1}} + \boldsymbol{e_{1}e_{2}} + \boldsymbol{e_{1}e_{3}} + \boldsymbol{e_{2}e_{1}} + \boldsymbol{e_{2}e_{2}} + \boldsymbol{e_{2}e_{3}} \\
&= 1 + \boldsymbol{e_{1}e_{2}} - \boldsymbol{e_{3}e_{1}} - \boldsymbol{e_{1}e_{2}} + 1 + \boldsymbol{e_{2}e_{3}} \\
&= 2 + \boldsymbol{e_{2}e_{3}} - \boldsymbol{e_{3}e_{1}}
\end{aligned}$$
and we can see that the scalar part, 2, is the regular dot product of $\boldsymbol{A}$ and $\boldsymbol{B}$. The bivector part $(\boldsymbol{e_{2}e_{3}} - \boldsymbol{e_{3}e_{1}})$, is actually the cross product of $\boldsymbol{A}$ and $\boldsymbol{B}$, expressed as a plane of rotation.

More info:

+ <https://en.wikipedia.org/wiki/Multivector>
+ <https://en.wikipedia.org/wiki/Clifford_algebra>
+ <https://people.kth.se/~dogge/clifford/files/clifford.pdf>
