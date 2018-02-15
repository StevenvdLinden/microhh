import numpy

# Get number of vertical levels and size from .ini file
with open('steven.ini') as f:
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
q_bg = numpy.zeros(numpy.size(z))

c0 = numpy.zeros(numpy.size(z))
c1 = numpy.zeros(numpy.size(z))
c2 = numpy.zeros(numpy.size(z))
pl0 = numpy.zeros(numpy.size(z))
pl1 = numpy.zeros(numpy.size(z))
fdn_above_0 = numpy.zeros(numpy.size(z))
fdn_above_1 = numpy.zeros(numpy.size(z))

u [:] = 0.
ug[:] = 0.

for k in range(rowl):
  if(z[k] <= 600.):
    th[k] = 300.
    q_bg[k] = 4.0e-3
  if(z[k] > 600.):
      #th[k] = 300. + dthetadz*(z[k]-100.)
    th[k] = 293. + dthetadz*(z[k]-100.)
    q_bg[k] = 4.0e-3

A = numpy.loadtxt('terms_Chiel')

for k in range(radmax):
  c0[k] = A[k][0]
  c1[k] = A[k][1]
  c2[k] = A[k][2]
  pl0[k] = A[k][3]
  pl1[k] = A[k][4]
  fdn_above_0[k] = A[k][5]
  fdn_above_1[k] = A[k][6]

# write the data to a file
proffile = open('steven.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s} {5:^20s} {6:^20s} {7:^20s} {8:^20s} {9:^20s} {10:^20s} {11:^20s} \n'.format('z','th','u','ug','q_bg','c0','c1','c2','pl0','pl1','fdn_above_0','fdn_above_1'))
for k in range(rowl):
  proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E} {6:1.14E} {7:1.14E} {8:1.14E} {9:1.14E} {10:1.14E} {11:1.14E} \n'.format(z[k], th[k], u[k], ug[k], q_bg[k], c0[k], c1[k], c2[k], pl0[k], pl1[k], fdn_above_0[k], fdn_above_1[k]))
proffile.close()

# write surface temperature
#timefile = open('steven.time','w')
#timefile.write('t     sbot[th] \n')
#timefile.write('0     265 \n')
#timefile.write('32400 265 \n')
#timefile.close()
