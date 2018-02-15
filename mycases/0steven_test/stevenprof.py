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

dthetadz = 0.

# set the height
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
th = numpy.zeros(numpy.size(z))
u = numpy.zeros(numpy.size(z))
ug = numpy.zeros(numpy.size(z))
#q_bg = numpy.zeros(numpy.size(z))

#c0 = numpy.zeros(numpy.size(z))
#c1 = numpy.zeros(numpy.size(z))
#c2 = numpy.zeros(numpy.size(z))
#pl0 = numpy.zeros(numpy.size(z))
#pl1 = numpy.zeros(numpy.size(z))
#fdn_above_0 = numpy.zeros(numpy.size(z))
#fdn_above_1 = numpy.zeros(numpy.size(z))

u [:] = 0.
ug[:] = 0.

for k in range(kmax):
  if(z[k] <= 100.):
    th[k] = 265.
  if(z[k] > 100.):
    th[k] = 265. + dthetadz*(z[k]-100.)

#for k in range(kmax):
# q_bg[k] = 0.
# c0[k] = 0.
# c1[k] = 0.
# c2[k] = 0.
# pl0[k] = 0.
# pl1[k] = 0.
# fdn_above_0[k] = 0.
# fdn_above_1[k] = 0.

# write the data to a file
proffile = open('steven.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} \n'.format('z','th','u','ug'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} \n'.format(z[k], th[k], u[k], ug[k]))
proffile.close()

# write surface temperature
#timefile = open('steven.time','w')
#timefile.write('t     sbot[th] \n')
#timefile.write('0     265 \n')
#timefile.write('32400 265 \n')
#timefile.close()
