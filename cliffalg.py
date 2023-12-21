###########

from itertools import combinations
from math import cos, sin, factorial, pi, sqrt
import re

class Algebra():
    """A class for real Clifford algebras.
    
    Takes a dictionary with names for basis vectors as keys and signatures as values.
    
    The algebra will always be over the field of real numbers, represented by float.
    
    May also be called with specific strings to get the following algebras
        'plane'      -->    {x: 1, y: 1}
        'space'      -->    {x: 1, y: 1, z: 1}
        'spacetime'  -->    {t: 1, x: -1, y: -1, z: -1}
        'complex'    -->    {i: -1}

    """

    def __init__(self, signature):
        """Constructs basis blades and their multiplication table."""
        
        if type(signature) is str:  # Translate named algebras into its dictionary.
            alg_names = {
                'plane':        {'x': 1, 'y': 1},
                'space':        {'x': 1, 'y': 1, 'z': 1},
                'spacetime':    {'t': 1, 'x': -1, 'y': -1, 'z': -1},
                'complex':      {'i': -1},
            }
            self._sign_dict = alg_names[signature]
        else:
            self._sign_dict = signature


        # Create a list of all basis blades from the names of the basis vectors.
        blade_list = []
        grade_list = [] # Keeps track of the grade of each basis blade.
        for n in range(len(signature)+1):
            blades = [list(comb) for comb in combinations(self._sign_dict.keys(), n,)]
            blade_list += blades
            grade_list += [n]*len(blades)
        
        # Construct multiplication table
        mult_table = []
        for b1 in blade_list:
            row = []
            for b2 in blade_list:
                commute_factor = 1
                blade = b1 + b2     # Concatenate the blade factors.

                # Bubble sort the vectors until the blade becomes a basis blade.
                i = 0
                while blade not in blade_list:
                    if i > len(blade)-2:  # Start over at 0 when the end of the list is reached.
                        i = 0

                    # Switch neighbors if the earlier one has an 'alphabetically greater' name.
                    if blade[i] > blade[i+1]:   # On each transposition, switch the sign of the coefficient.
                        commute_factor *= -1             # This is the anticommutative property of the Clifford product.
                        blade[i], blade[i+1] = blade[i+1], blade[i]
                    
                    # Two of the same vectors directly beside each other yield the signature.
                    elif blade[i] == blade[i+1]:
                        commute_factor *= self._sign_dict[blade.pop(i)]  # The first copy is removed here.
                        blade.pop(i)                            # The second copy is removed here.
                    
                    i += 1  # Step forward in the list of vectors. 
                
                # Only store the coefficient and the index of the basis blade in the table.
                row.append((commute_factor, blade_list.index(blade)))    
            mult_table.append(tuple(row))

        # Construct instance variables as tuples for extra stability.
        self._blade_names = tuple(''.join(lis) for lis in blade_list)
        self._grades = tuple(grade_list)
        self._mult_table = tuple(mult_table)
    
    def sign_dict(self):
        """Returns the dictionary of basis blades and their signatures."""
        return self._sign_dict

    def __repr__(self):
        """Returns a string from the dict of basis vectors and their signatures."""

        return f'Algebra({self._sign_dict})'
    
    def __eq__(self, other):
        """Two algebras are equal if and only if they have the same set of 
        basis vectors and the respective signatures are equal."""

        return self.sign_dict() == other.sign_dict()
    
    def blade_names(self):
        """Returns a tuple of strings with the names of all basis blades.
        
        Worst case time complexity: O(b) = O(2^v), where v is the number of basis vectors in
        the algebra and b = 2^v is the number of basis blades in the algebra.
        """
        return self._blade_names

    def blades(self):
        """Returns a tuple of all basis blades in this algebra.
        
        All blades are given in the form of multivectors where one component is 1
        and all the others are 0. E.g. the first element of the tuple will be 
        the object Multivector(1, self), not int(1).

        Worst case time complexity: O(b^2) = O(4^v), since blade_names are also
        iterated over when creating a new multivector.
        """
        return tuple(Multivector({name: 1}, self) for name in self.blade_names())


class Multivector():
    """A class for multivectors in a given Clifford algebra.
    
    Can be constructed from:

        (#) A list or tuple of components. The basis blades are ordered by grade and
            alphabetically within each grade. E.g. 1, x, y, z, xy, xz, yz, xyz.

        (#) A dictionary with names of basis blades as keys and components as values.
            Unspecified components are assumed to be 0.
            
        (#) An int of float ('real scalar'). This results in a scalar in the algebra.
            
            E.g. all the following constructs gives a zero in an algebra called alg:
            
                Multivector({'': 0}, alg) == Multivector({}, alg) == Multivector(0, alg)
    
    
    Arithmetic operations:

        (#) Supports addition (+), subtraction (-) and Clifford product (*) with other 
            multivectors of the same algebra or with real scalars i.e. ints or floats.
            
        (#) A multivector may be divided (/) by a real scalar but no type may be 
            divided by a multivector.

        (#) Negation is supported.

    """

    def __init__(self, components, algebra):
        """Creates a new multivector object.
        
        Worst case time complexities:
        
            (#) list or tuple: O(b),

            (#) dict: O(b + k*d),

            (#) int or float: O(b + d),

        where b is the number of blades in the algebra, k is the number of 
        keys in the dictionary and d is the number of digits per component. 
        """
        self._algebra = algebra

        # Create component tuple based on input or raise errors if the input is bad.
        if (type(components) in (tuple, list)) and all(type(comp) in (int, float) for comp in components):
            if len(components) == len(algebra.blade_names()):
                self._comp_tup = tuple(components)
            else:
                raise IndexError
        else:
            comp_list = [0]*len(algebra.blade_names())
        
            if (type(components) is dict) and all(type(comp) in (int, float) for comp in components.values()):
                if all(blade in algebra.blade_names() for blade in components.keys()):
                    for blade, comp in components.items():
                        comp_list[algebra.blade_names().index(blade)] = comp
                else:
                    raise KeyError
            elif type(components) in (int, float):
                comp_list[0] = components
            else:
                raise TypeError

            self._comp_tup = tuple(comp_list)
    
    def algebra(self):
        """Returns the algebra object that this belongs to.

        Worst case time complexity: O(1).
        """
        return self._algebra

    def comps(self, *grades):
        """Return the components of this multivector as a tuple.
        
        Which grades to include can be specified. Default is all.

        Worst case time complexity: O(b) with specified grades, otherwise O(1).
        """
        if len(grades) == 0:
            return self._comp_tup
        else:
            return tuple(comp for comp, grade in zip(self._comp_tup, self.algebra()._grades) if grade in grades)
            

    def __repr__(self):
        """Returns a str of this as a linear combination of basis blades."""

        indices = [i for i, comp in enumerate(self.comps()) if comp != 0]  # Only show nonzero components
        if len(indices) == 0:
            return '0'
        else:
            str1 = ' + '.join([f'{self.comps()[i]}{self.algebra().blade_names()[i]}' for i in indices])
            str2 = re.sub(r'\+ -','- ', str1)            # Trim plus-signs in front of negative components.
            str3 = re.sub(r'(?<![^ -])1(?=[^0-9 .])', '', str2)    # Trim redundant 1-coefficients.
            return str3
    
    def __eq__(self, other):
        """Two multivectors are equal if and only if they belong to the same
        algebra and their components are equal.
        
        Worst case time complexity: O(b).
        """

        return self.comps() == other.comps()
    
    
    def __neg__(self):
        """Negates this multivector by multiplying all components by -1.
        
        Worst case time complexity: O(b).
        """
        
        components = tuple(-comp for comp in self.comps())
        return Multivector(components, self.algebra())

    def __add__(self, other):
        """Returns the sum of this multivector with another object.
        
        The object may be a multivector from the same algebra or a real scalar.
        This method is automatically called when evaluating (self + other).

        Worst case time complexity: O(b*d).
        """
        if type(other) is not Multivector:              # Convert real scalars to multivector 
            other = Multivector(other, self.algebra())
        components = tuple(self.comps()[i] + other.comps()[i] for i in range(len(self.comps())))
        return Multivector(components, self.algebra())
    
    def __radd__(self, other):
        """Returns the sum of this multivector with another object.
        
        The object may be a multivector from the same algebra or a real scalar.
        This method is automatically called when evaluating (other + self).

        Worst case time complexity: O(b*d).
        """
        if type(other) is not Multivector:              # Convert real scalars to multivector 
            other = Multivector(other, self.algebra())
        components = tuple(self.comps()[i] + other.comps()[i] for i in range(len(self.comps())))
        return Multivector(components, self.algebra())

    def __sub__(self, other):
        """Returns the arithmetic difference between this multivector and another object.
        
        The object may be a multivector from the same algebra or a real scalar.
        This method is automatically called when evaluating (self - other).

        Worst case time complexity: O(b*d).
        """
        if type(other) is not Multivector:              # Convert real scalars to multivector 
            other = Multivector(other, self.algebra())
        components = tuple(self.comps()[i] - other.comps()[i] for i in range(len(self.comps())))
        return Multivector(components, self.algebra())
    
    def __rsub__(self, other):
        """Returns the arithmetic difference between another object and this multivector.
        
        The object may be a multivector from the same algebra or a real scalar.
        This method is automatically called when evaluating (other - self).

        Worst case time complexity: O(b*d).
        """
        if type(other) is not Multivector:              # Convert real scalars to multivector 
            other = Multivector(other, self.algebra())
        components = tuple(other.comps()[i] - self.comps()[i] for i in range(len(self.comps())))
        return Multivector(components, self.algebra())
    
    def __mul__(self, other):
        """Returns the product of this with another object, with this on the left. 
        
        The object may be a multivector from the same algebra or a real scalar.
        This method is automatically called when evaluating (self*other).

        Worst case time complexity: O(b^2*d^2) when other is a multivector.
        O(b*d^2) when other is a real scalar.
        """
        if type(other) in (int, float):         # Multiply components directly if other is a real scalar.
            components = tuple(other*comp for comp in self.comps())
            return Multivector(components, self.algebra())
        components = [0]*len(self.comps())    # List to store terms
        for i, comp1 in enumerate(self.comps()):
            if comp1 != 0:                      # Skip if one of the factors is zero
                for j, comp2 in enumerate(other.comps()):
                    if comp2 != 0:              # Skip if one of the factors is zero
                        commute_factor, blade_index = self.algebra()._mult_table[i][j]   # Extract product of blade i with blade j.
                        components[blade_index] += commute_factor*comp1*comp2            # Add the term to the correct index.

        return Multivector(components, self.algebra())

    
    def __rmul__(self, other):
        """Returns the product of this with another object, with this on the right. 
        
        The object may be a multivector from the same algebra or a real scalar.
        This method is automatically called when evaluating (other*self).

        Worst case time complexity: O(b^2*d^2) when other is a multivector.
        O(b*d^2) when other is a real scalar.
        """
        if type(other) in (int, float):         # Multiply components directly if other is a real scalar.     
            components = tuple(other*comp for comp in self.comps())
            return Multivector(components, self.algebra())

        components = [0]*len(self.comps())      # List to store terms
        for i, comp1 in enumerate(other.comps()):
            if comp1 != 0:                      # Skip if one of the factors is zero
                for j, comp2 in enumerate(self.comps()):
                    if comp2 != 0:              # Skip if one of the factors is zero
                        commute_factor, blade_index = self.algebra()._mult_table[i][j]   # Extract product of blade i with blade j.
                        components[blade_index] += commute_factor*comp1*comp2            # Add the term to the correct index.

        return Multivector(components, self.algebra())

    def __truediv__(self, x):
        """Divide this multivector by a real scalar x.
        
        Worst case time complexity: O(b*d^2).
        """

        components = tuple(comp/x for comp in self.comps())
        return Multivector(components, self.algebra())
    

    def __pow__(self, p):
        """The nth power of this multivector, where p is a non negative integer.

        Worst case time complexity: O(b^2*d^2*p(p-1)/2).

        The triangluar p-factor comes from the fact that the number of digits in
        the components of the subproduct increases by d each step.
        """
        
        sub_prod = Multivector(1, self.algebra())
        for _ in range(p):
            sub_prod = sub_prod*self
        return sub_prod
    
    def exp(self, taylor_terms=16, _disable_trig=False):
        """Returns a new multivector which is the exponential (base e) of this.
        
        Calculated with the trig identity exp(ix) = cos(x) + isin(x) if this multivector
        squares to -1.
        
        Otherwise calculated with the taylor series exp(self) = 1 + self + (self**2)/2 +...
        The number of terms used may be specified. Default is 16.

        Worst case time complexities:
        
            (#) trig identity: O(b^2*d^2),

            (#) taylor series: O(b^2*d^2*t^2),

        where t=taylor_terms
        """

        # Trig identity
        if not _disable_trig:  # Can be disabled for debugging.
            square = self**2                # Check if this squares to a negative real scalar.
            sq_scal = square.comps(0)[0]
            if sq_scal < 0 and all(comp == 0 for comp in square.comps()[1:]):
                magn = sqrt(-sq_scal)       # Normalize and use normalization factor as angle.
                return (cos(magn) + self/magn*sin(magn))

        # Taylor series
        result = 1
        for n in range(1, taylor_terms):
            result += (self**n)/factorial(n)
        return result
            
    
    def grade_proj(self, k):
        """Returns a new multivector that only keeps the componets of grade k.
        
        Worst case time complexity: O(b).
        """
        components = tuple(comp if grade == k else 0 for comp, grade in zip(self.comps(), self.algebra()._grades))
        return Multivector(components, self.algebra())
    

def test():
    R3 = Algebra('space')
    x, y, z, xy, xz, yz, xyz = R3.blades()[1:]

    u = Multivector(-1,R3)
    v = 2.1*xy
    w = 1 + 2.21*xy
    print(u)
    print(v)
    print(w)

    A = Multivector([-1., -2., 3., -5., 7., 11., 13., -4.], R3)
    B = Multivector([-7.,1.,16.,-24.,93.,-1.2,2.2,31.], R3)
    print(A*B)

if __name__ == '__main__':
    test()
