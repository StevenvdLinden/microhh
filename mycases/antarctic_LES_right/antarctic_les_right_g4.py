import numpy

# Get number of vertical levels and size from .ini file
with open('antarctic_test_new_g4coeff.ini') as f:
  for line in f:
    if(line.split('=')[0]=='ktot'):
      kmax = int(line.split('=')[1])
    if(line.split('=')[0]=='zsize'):
      zsize = float(line.split('=')[1])

# set the height
dz = zsize / kmax

# LELIJKE METHODE, MAAR NA BEPALING dz MAAK ARRAY 'KUNSTMATIG' LANGER. DEZE WAARDES WORDEN LATER TOCH NIET INGELEZEN

radmax = 55
# kmax = radmax

dthetadz = 0.

rowl = numpy.max([kmax,radmax])

# set the height
z = numpy.linspace(0.5*dz, dz*rowl-0.5*dz, rowl) #deze is ook aangepast, gebruik niet de oude zsize!
th = numpy.zeros(numpy.size(z))
u = numpy.zeros(numpy.size(z))
ug = numpy.zeros(numpy.size(z))
wls = numpy.zeros(numpy.size(z))
q_bg = numpy.zeros(numpy.size(z))

c0 = numpy.zeros(numpy.size(z))
c1 = numpy.zeros(numpy.size(z))
c2 = numpy.zeros(numpy.size(z))
pl0 = numpy.zeros(numpy.size(z))
pl1 = numpy.zeros(numpy.size(z))
fdn_above_0 = numpy.zeros(numpy.size(z))
fdn_above_1 = numpy.zeros(numpy.size(z))

u [:] = 15.
ug[:] = 15.

for k in range(rowl):
  if(z[k] <= 90.):
    th[k] = 220 + (25. / 90.) * z[k]
    q_bg[k] = 3.7e-4
  if(z[k] > 90.):
      #th[k] = 300. + dthetadz*(z[k]-100.)
    th[k] = 245. + dthetadz*(z[k]-100.)
    q_bg[k] = 3.7e-4

for k in range(rowl):
    wls[k] = -0.01 * (z[k] / 100.)

A = numpy.loadtxt('terms_G4_top100m_microhh')

for k in range(radmax):
  c0[k] = A[k][0]
  c1[k] = A[k][1]
  c2[k] = A[k][2]
  pl0[k] = A[k][3]
  pl1[k] = A[k][4]
  fdn_above_0[k] = A[k][5]
  fdn_above_1[k] = A[k][6]

# write the data to a file
proffile = open('antarctic_test_new_g4coeff.prof','w')
proffile.write('{0:^16s} {1:^16s} {2:^16s} {3:^16s} {4:^16s} {5:^16s} {6:^16s} {7:^16s} {8:^13s} {9:^13s} {10:^13s} {11:^13s} {12:^13s}\n'.format('z','th','u','ug','wls','q_bg','c0','c1','c2','pl0','pl1','fdn_above_0','fdn_above_1'))
for k in range(rowl):
  proffile.write('{0:1.10E} {1:1.10E} {2:1.10E} {3:1.10E} {4:1.10E} {5:1.10E} {6:1.7E} {7:1.7E} {8:1.7E} {9:1.7E} {10:1.7E} {11:1.7E} {12:1.7E}\n'.format(z[k], th[k], u[k], ug[k], wls[k], q_bg[k], c0[k], c1[k], c2[k], pl0[k], pl1[k], fdn_above_0[k], fdn_above_1[k]))
proffile.close()

# write surface temperature
#timefile = open('steven.time','w')
#timefile.write('t     sbot[th] \n')
#timefile.write('0     265 \n')
#timefile.write('32400 265 \n')
#timefile.close()
