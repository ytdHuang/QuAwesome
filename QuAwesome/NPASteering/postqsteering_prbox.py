import qutip as qt
import numpy as np
import matplotlib.pyplot as plt
import steering as st
import picos as pic


# PR box correlations
def prbox(a,b,x,y):
    if (a+b)%2 == x*y:
        return 1/2
    else:
        return 0

prb = {}
for a in range(2):
    for b in range(2):
        for x in range(2):
            for y in range(2):
                prb[str(a)+str(b)+str(x)+str(y)] = prbox(a,b,x,y)


# define assemblage
assemb = {}
for a in range(2):
    for b in range(2):
        for x in range(2):
            for y in range(2):
                assemb[str(a)+str(b)+str(x)+str(y)] = prb[str(a)+str(b)+str(x)+str(y)]*qt.qeye(2)/2

x=1;y=1
rho_r = 0
for a in range(2):
    for b in range(2):
        rho_r += assemb[str(a)+str(b)+str(x)+str(y)]
rho_r

assemb_a = {}
for x in range(2):
    for a in range(2):
        q=0
        for b in range(2):
            q+=assemb[str(a)+str(b)+str(x)+str(y)]
        assemb_a[str(a)+str(x)]=q
        
assemb_b = {}
for y in range(2):
    for b in range(2):
        q=0
        for a in range(2):
            q+=assemb[str(a)+str(b)+str(x)+str(y)]
        assemb_b[str(b)+str(y)]=q

# assemblage as picos object
r = pic.Constant('r',rho_r.full())

s={}
for x in range(2):
    for y in range(2):
        for a in range(2):
            for b in range(2):
                s[str(a)+str(b)+str(x)+str(y)]=pic.Constant('s_'+str(a)+str(b)+str(x)+str(y),
                                                                 assemb[str(a)+str(b)+str(x)+str(y)].full())
sa={}
for x in range(2):
    for a in range(2):
        sa[str(a)+str(x)]=pic.Constant('sa_'+str(a)+str(x), assemb_a[str(a)+str(x)].full())

sb={}
for y in range(2):
    for b in range(2):
        sb[str(b)+str(y)]=pic.Constant('sb_'+str(b)+str(y), assemb_b[str(b)+str(y)].full())

# defining varibles X_i

X = [pic.HermitianVariable('X{0}'.format(i) , (2,2)) for i in range(8)]

# constructing chi matrix

r1 = (r & sa['00'] & sa['01'] & 
      sb['00'] & s['0000'] & s['0010'] &
      sb['01'] & s['0001'] & s['0011']
     )

r2 = (sa['00'] & sa['00'] & X[0] & 
      s['0000'] & s['0000'] & X[1] &
      s['0001'] & s['0001'] & X[2]
     )


r3 = (sa['01'] & X[0] & sa['01'] & 
      s['0010'] & X[1] & s['0010'] &
      s['0011'] & X[2] & s['0011']
     )


r4 = (sb['00'] & s['0000'] & s['0010'] & 
      sb['00'] & s['0000'] & s['0010'] &
      X[3] & X[4] & X[5]
     )


r5 = (s['0000'] & s['0000'] & X[1] & 
      s['0000'] & s['0000'] & X[1] &
      X[4] & X[4] & X[6]
     )


r6 = (s['0010'] & X[1] & s['0010'] & 
      s['0010'] & X[1] & s['0010'] &
      X[5] & X[7] & X[5]
     )


r7 = (sb['01'] & s['0001'] & s['0011'] & 
      X[3] & X[4] & X[5] &
      sb['01'] & s['0001'] & s['0011']
     )


r8 = (s['0001'] & s['0001'] & X[2] & 
      X[4] & X[4] & X[7] &
      s['0001'] & s['0001'] & X[2]
     )


r9 = (s['0011'] & X[2] & s['0011'] & 
      X[5] & X[6] & X[5] &
      s['0011'] & X[2] & s['0011']
     )


chi = r1//r2//r3//r4//r5//r6//r7//r8//r9

post = pic.Problem()
post.add_constraint(chi>>0)

post.set_objective('find')
res = post.solve(solver = 'mosek',verbosity=1)

