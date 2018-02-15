import numpy
from pylab import *

# set the height
kmax  = 96

gf = 1.072
nu = 0.002
G  = 5.
f  = 0.04

ust = G/18
vsc = nu/ust
dz0 = vsc/4
gamma = 1/sqrt(2*nu/f)

z = zeros(kmax)
u = zeros(kmax)
v = zeros(kmax)
s = zeros(kmax)
ug = zeros(kmax)
vg = zeros(kmax)

ug[:] = G
vg[:] = 0.

#u[:] = ug[:]
#v[:] = vg[:]

z[0] = dz0

for k in range(kmax-1):
    if(k==0):
        z[k+1] = z[k] + gf*(z[k])
    else:
        z[k+1] = z[k] + gf*(z[k]-z[k-1])

# analytical solution as the starting profile to reduce run time
for k in range(kmax):
  u[k] = ug[k]*(1. - exp(-gamma*z[k]) * cos(gamma*z[k]))
  v[k] = ug[k]*(     exp(-gamma*z[k]) * sin(gamma*z[k]))

# write the data to a file
proffile = open('ekman.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s}\n'.format('z','u','v','ug','vg','s'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E}\n'.format(z[k], u[k], v[k], ug[k], vg[k], s[k]))
proffile.close()

