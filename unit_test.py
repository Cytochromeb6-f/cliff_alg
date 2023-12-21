from cliffalg import Algebra, Multivector

from random import choice, random


# Used to create multivectors with (uniformly) random float coefficients in a specified interval.
def randmv(a, b, alg):
    return Multivector({name: a + (b-a)*random() for name in alg.blade_names()}, alg)

# Test named algebras.
def named_algebras():
    algebra_dict = {
        'plane': {'x': 1, 'y': 1},
        'space': {'x': 1, 'y': 1, 'z': 1},
        'spacetime': {'t': 1, 'x': -1, 'y': -1, 'z': -1},
        'complex': {'i': -1}
    }

    for name, sg_dict in algebra_dict.items():
        assert Algebra(name)                == Algebra(sg_dict)
        assert Algebra(name).sign_dict()    == sg_dict
        assert str(Algebra(name))           == 'Algebra(' + str(sg_dict) + ')'


# Test creating multivectors.
def creating_multivectors(repeat):
    for alg in (Algebra('plane'), Algebra('space'), Algebra('spacetime'), Algebra('complex')):
        for _ in range(int(repeat)):
            N = len(alg.blade_names())
            a = 20*random()-10          # Random float between -10 and 10
            assert Multivector(a, alg) == Multivector({'': a}, alg) == Multivector([a]+[0]*(N-1), alg)


# Test projection onto a specified grade.
def grade_projection(repeat):
    for alg in (Algebra('plane'), Algebra('space'), Algebra('spacetime'), Algebra('complex')):
        for _ in range(int(repeat)):
            u = randmv(-10, 10, alg)
            
            for k in range(len(alg.sign_dict())+1):
                comp_list = [comp if len(name) == k else 0 for comp, name in zip(u.comps(), alg.blade_names())]    
                assert u.grade_proj(k) == Multivector(comp_list, alg)


# Test getting the basis blades of an algebra.
def blades():
    blade_dict = {
        'plane': ('', 'x', 'y', 'xy'),
        'space': ('', 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'),
        'spacetime': ('', 't', 'x', 'y', 'z', 'tx', 'ty', 'tz', 'xy', 'xz', 'yz', 'txy', 'txz', 'tyz', 'xyz','txyz'),
        'complex': ('', 'i')
    }

    for alg_name in blade_dict.keys():
        alg = Algebra(alg_name)
        assert alg.blade_names() == blade_dict[alg_name]

        for blade, name in zip(alg.blades(), blade_dict[alg_name]):
            assert blade == Multivector({name: 1}, alg)


# Test string representation of multivectors.
def string_representation(repeat):
    # A multivector is represented by a linear combination of coefficient-bladename pairs.

    for alg in (Algebra('plane'),):# Algebra('space'), Algebra('spacetime'), Algebra('complex')):
        for _ in range(int(repeat)):
        
            # Construct a multivector with random coefficients that are either -2,-1,0,1,2 or a random float.
            comp_list = []
            for _ in range(len(alg.blade_names())):
                alternatives = (-2, -1, 0, 1, 2, 20*random()-10)
                comp_list.append(choice(alternatives))
            u = Multivector(comp_list, alg)
            
            # Construct a string of what the multivector should look like.
            if all(coeff == 0 for coeff in u.comps()):
                correct = '0'               # Should return '0' if all coefficients are zero.
            else:
                correct = ''
                for coeff, b_name in zip(u.comps(), alg.blade_names()):
                    
                    # Skip whole term if coefficient is zero.
                    if coeff == 0:          
                        continue
                    
                    # A scalar term is always just its coefficient.
                    if b_name == '':
                        correct += str(coeff)
                        continue
                    
                    # First nonzero, nonscalar term.
                    if correct == '':
                        if coeff > 0:       # Operator is '' for positive coefficients.
                            op = ''
                        else:               # Operator is '-' for negative coefficients.
                            op = '-'
                    
                    # Subsequent terms.
                    else:
                        if coeff > 0:       # Operator is ' + ' for positive coefficients.
                            op = ' + '
                        else:               # Operator is ' - ' for negative coefficients.
                            op = ' - '
                    
                    if abs(coeff) == 1: # Skip redundant 1-coefficients.
                        correct += f'{op}{b_name}'
                    else:
                        correct += f'{op}{abs(coeff)}{b_name}'

            assert str(u) == correct   # Compare the strings.


# Test arithmetic.
def arithmetic(repeat, TOL):
    for alg in (Algebra('plane'), Algebra('space'), Algebra('spacetime'), Algebra('complex')):    
        for _ in range(int(repeat)):
            
            # Test with u as a real scalar and u as a multivector.
            for u in (20*random()-10, randmv(-10, 10, alg)):
                v = randmv(-10, 10, alg)
                w = randmv(-10, 10, alg)
                c = 20*random()-10

                # Basic tests.
                assert u == u       # __eq__
                assert u+v == v+u   # __add__, __radd__
                assert u-v==-(v-u)  # __neg__, __sub__, __rsub__


                # Multiplication tests. Small numerical errors are tolerated.
                diff = [            # Clifford product should satisfy the axioms for a bilinear form.
                    (u+v)*w - (u*w+v*w),    # case 0    Right distributivity  __mul__
                    u*(v+w) - (u*v+u*w),    # case 1    Left distributivity   __rmul__
                    (c*u)*v - c*(u*v),      # case 2    Scalar compatability
                    u*(c*v) - c*(u*v),      # case 3    
                                    # Extra tests
                    -u*(v+w) + (u*v+u*w),   # case 4    __neg__
                    (u-v)*w - (u*w-v*w),    # case 5    __sub__, __rsub__
                    u*(v/c) - (u*v)/c,      # case 6    __truediv__
                    (1/c)*v - v/c           # case 7    
                ]
                for i, mv in enumerate(diff):
                    
                    # Make sure that the new multivectors stay in the same algebra. 
                    assert mv.algebra() == alg      # All Clifford algebras are closed under arithmetic.

                    # The maximum deviation should be less than the desired tolerance.
                    devs = [abs(comp) for comp in mv.comps()]
                    try:
                        assert max(devs) < TOL
                    except AssertionError:
                        print(f'arithmetic | TOL: {TOL}, test case: {i}, dev: {max(devs)}, {alg}')


# Test powers
def powers(max_power, repeat, TOL):
    for alg in (Algebra('plane'), Algebra('space'), Algebra('spacetime'), Algebra('complex')):
        for _ in range(int(repeat)):
            u1 = randmv(-10, 10, alg)

            u = Multivector(1, alg)
            for p in range(max_power+1):    
                assert u == u1**p

                u = u*u1        # Left multiplication. This is how the method is implemented.
            
            u = Multivector(1, alg)
            for p in range(max_power+1):        
                try:
                    assert u == u1**p
                except AssertionError:
                    dev = []
                    for comp, diff_comp in zip(u.comps(), (u-u1**p).comps()):
                        if comp == 0:                   # Use absolute deviation if numerator is zero.
                            dev.append(abs(diff_comp))
                        else:                           # Otherwise use relative deviation.
                            dev.append(abs(diff_comp/comp))
                        
                    if max(dev) > TOL:
                        print(f'    powers | TOL: {TOL}, power: p = {p}, dev: {max(dev)}, {alg}')

                u = u1*u        # Right multiplication. Should be the same algebraically.
                    



if __name__ == '__main__':

    named_algebras()

    creating_multivectors(repeat=50)

    grade_projection(repeat=50)

    blades()
    
    string_representation(repeat=1e3)



    rep = 1e2       # Quick
    #rep = 1e3       # Exhaustive

    tol = 1e-9      # Safe
    #tol = 1e-12     # Too intolerant

    arithmetic(repeat=rep, TOL=tol)
    powers(max_power=9, repeat=rep, TOL=tol)

   
