from sympy.solvers import solve
from sympy import Symbol

eeff = Symbol('eeff')
f = Symbol('f')
e0 = Symbol('e0')
e1 = Symbol('e1')
v = Symbol('v')

print( solve((eeff-e0)/(eeff+2.0*e0+v*(eeff-e0))-f*(e1-e0)/(e1+2.0*e0+v*(eeff-e0)),eeff) )
