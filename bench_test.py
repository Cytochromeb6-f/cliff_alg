
from cliffalg import Algebra, Multivector

from time import time
from random import random

from unit_test import randmv

# Investigates precision losses when multiplying multivectors.
def product_precision(max_power):
    print(f'\n  Product precision test')

    u = randmv(-10, 10, Algebra('space'))
    
    prod = [1, 1] 
    for p in range(max_power+1):
        prod[0] = u*prod[0]         # Left multiplication
        prod[1] = prod[1]*u         # Right multiplication
        dev_list = [abs(comp) for comp in (prod[0]-prod[1]).comps()]
        abs_dev = max(dev_list)
        rel_dev = max([dev/comp for dev, comp in zip(dev_list, prod[0].comps())])

        print(f'p = {p:{2}}, abs_dev = {abs_dev:{8}.{3}}, rel_dev = {rel_dev:{8}.{3}}')
    


# Compares speed of creating multivectors with dict, list and tuple.
def creation_speed(n):
    print(f'\n  Creation comparison, n = {n}')
    for alg in (Algebra('plane'), Algebra('space'), Algebra('spacetime')):
        print(alg)
        dicts = []
        tuples = []
        lists = []
        for _ in range(int(n)):
            di = {name: 20*random()-10 for name in alg.blade_names()}
            li = list(di.values())
            tu = tuple(li)
            dicts.append(di)
            lists.append(li)
            tuples.append(tu)

        t1 = time()
        for tu in tuples:
            Multivector(tu,alg)
        t2 = time()
        print(f'Tuples : {t2-t1} s')

        t1 = time()
        for li in lists:
            Multivector(li, alg)
        t2 = time()
        print(f'Lists  : {t2-t1} s')

        t1 = time()
        for di in dicts:
            Multivector(di, alg)
        t2 = time()
        print(f'Dicts  : {t2-t1} s')

# Compares trig and taylor exponentiation
def exp_speed(n):
    print(f'\n  Exp speed comparison, n = {n}')
    R3 = Algebra('space')
    bivectors = []          # Test on randomly generated bivectors.
    for _ in range(int(n)):
        a, b, c = random(), random(), random()
        bivectors.append(Multivector((0, 0, 0, 0, a, b, c, 0), R3))
    t1 = time()
    for b in bivectors:
        b.exp()
    t2 = time()
    print(f'Trig  : {t2-t1} s')
    for terms in (2,4,8,12,16):
        t1 = time()
        for b in bivectors:
            b.exp(taylor_terms=terms, _disable_trig=True)
        t2 = time()
        print(f'Taylor: {t2-t1} s  (using {terms} terms)')

if __name__ == '__main__':

    product_precision(20)

    creation_speed(5e3)

    exp_speed(500)
