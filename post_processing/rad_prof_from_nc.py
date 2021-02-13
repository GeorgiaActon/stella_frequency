

from scipy.io import netcdf
import numpy as np
import numpy.matlib


tave = 1000

basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/2nd_deriv/T1/'
basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/2nd_deriv/rho_scan2/r0.001/'
basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/fg_drive/0.001hr/'
right_file  = basedir + 'right.out.nc'
center_file = basedir + 'center.out.nc'
left_file   = basedir + 'left.out.nc'


right_nc  = netcdf.netcdf_file(right_file,'r')
center_nc = netcdf.netcdf_file(center_file,'r')
left_nc   = netcdf.netcdf_file(left_file,'r')


def read_stella_float(infile, var):

  import numpy as np

  try:
    #print('a')
    #arr = np.copy(infile.variables[var][:])
    arr = infile.variables[var][:]
    #print('b')
    flag = True
  except KeyError:
    print('INFO: '+var+' not found in netcdf file')
    arr =np.arange(1,dtype=float)
    flag = FLAG

  return arr, flag

def phi_vs_t_to_x(infile,var,ny,nx):
# t ntube z kx ky ri

  avt, present = read_stella_float(infile,var)
  #print('c')
  avt_kxky = ny*nx*(avt[:,0,:,:,:,0] + 1j*avt[:,0,:,:,:,1])
  #print('d')
  arr = np.fft.ifft(avt_kxky,axis=2)
  #print('e')

  return arr

def mom_vs_t_to_x(infile,var,ny,nx):
  #in:  t nspec ntube z kx ky ri
  #out: t z kx ky

  avt, present = read_stella_float(infile,var)
  avt_kxky = ny*nx*(avt[:,0,0,:,:,:,0] + 1j*avt[:,0,0,:,:,:,1])
  arr = np.fft.ifft(avt_kxky,axis=2)

  return arr

print('0')
naky  = center_nc.dimensions['ky']
nakxl =   left_nc.dimensions['kx']
nakxc = center_nc.dimensions['kx']
nakxr =  right_nc.dimensions['kx']
ky  = np.copy(center_nc.variables['ky'][:])
kxc = np.copy(center_nc.variables['kx'][:])

t  = np.copy(center_nc.variables['t'][:])
nt = t.size

Lxc = 2.*np.pi/kxc[1]
dxc = Lxc / nakxc

zed  = np.copy(center_nc.variables['zed'][:])
nzed = zed.size
omp = ((nzed+1)/2) - 1
delzed = zed[1]-zed[0]

fac = 2*np.ones(naky)
fac[0] = 1

jacobl  = np.copy(  left_nc.variables['jacob'][:])
jacobc  = np.copy(center_nc.variables['jacob'][:])
jacobr  = np.copy( right_nc.variables['jacob'][:])

print('1')

dl_over_bl = np.squeeze(delzed*jacobl/sum(delzed*jacobl))
dl_over_bc = np.squeeze(delzed*jacobc/sum(delzed*jacobc))
dl_over_br = np.squeeze(delzed*jacobr/sum(delzed*jacobr))

dobl = np.transpose(np.matlib.tile(dl_over_bl,(naky,nakxl,1)))
dobc = np.transpose(np.matlib.tile(dl_over_bc,(naky,nakxc,1)))
dobr = np.transpose(np.matlib.tile(dl_over_br,(naky,nakxr,1)))

print('2')

# t spec x
pfluxl  = np.copy(  left_nc.variables['pflux_x'][:])
pfluxc  = np.copy(center_nc.variables['pflux_x'][:])
pfluxr  = np.copy( right_nc.variables['pflux_x'][:])

vfluxl  = np.copy(  left_nc.variables['vflux_x'][:])
vfluxc  = np.copy(center_nc.variables['vflux_x'][:])
vfluxr  = np.copy( right_nc.variables['vflux_x'][:])

qfluxl  = np.copy(  left_nc.variables['qflux_x'][:])
qfluxc  = np.copy(center_nc.variables['qflux_x'][:])
qfluxr  = np.copy( right_nc.variables['qflux_x'][:])

cout = open(basedir + 'left.fluxes_t','w')
cout.write('[1] t     ')
cout.write('[2] x     ')
cout.write('[3] flux_d')
cout.write('[4] flux_u')
cout.write('[5] flux_t')
cout.write('\n')
print('3')
for i in range (0, nt):
  for j in range (0, nakxl):
    cout.write('%e ' % t[i])
    cout.write('%e ' % (dxc*j))
    cout.write('%e ' % pfluxl[i,0,j])
    cout.write('%e ' % vfluxl[i,0,j])
    cout.write('%e ' % qfluxl[i,0,j])
    cout.write('\n')
  cout.write('\n')
cout.close()

cout = open(basedir + 'center.fluxes_t','w')
cout.write('[1] t     ')
cout.write('[2] x     ')
cout.write('[3] flux_d')
cout.write('[4] flux_u')
cout.write('[5] flux_t')
cout.write('\n')
print('4')
for i in range (0, nt):
  for j in range (0, nakxc):
    cout.write('%e ' % t[i])
    cout.write('%e ' % (dxc*j))
    cout.write('%e ' % pfluxc[i,0,j])
    cout.write('%e ' % vfluxc[i,0,j])
    cout.write('%e ' % qfluxc[i,0,j])
    cout.write('\n')
  cout.write('\n')
cout.close()

cout = open(basedir + 'right.fluxes_t','w')
cout.write('[1] t     ')
cout.write('[2] x     ')
cout.write('[3] flux_d')
cout.write('[4] flux_u')
cout.write('[5] flux_t')
cout.write('\n')
print('5')
for i in range (0, nt):
  for j in range (0, nakxr):
    cout.write('%e ' % t[i])
    cout.write('%e ' % (dxc*j))
    cout.write('%e ' % pfluxr[i,0,j])
    cout.write('%e ' % vfluxr[i,0,j])
    cout.write('%e ' % qfluxr[i,0,j])
    cout.write('\n')
  cout.write('\n')
cout.close()

tind=nt-1
for i in range (0, nt):
  if(t[i]> tave):
    tind = i
    break
    
print(str(tind) + '  ' + str(nt))

print('6')

plave = np.mean(pfluxl[tind:nt,0,:],0)
vlave = np.mean(vfluxl[tind:nt,0,:],0)
qlave = np.mean(qfluxl[tind:nt,0,:],0)

pcave = np.mean(pfluxc[tind:nt,0,:],0)
vcave = np.mean(vfluxc[tind:nt,0,:],0)
qcave = np.mean(qfluxc[tind:nt,0,:],0)

prave = np.mean(pfluxr[tind:nt,0,:],0)
vrave = np.mean(vfluxr[tind:nt,0,:],0)
qrave = np.mean(qfluxr[tind:nt,0,:],0)

print('7')

cout = open(basedir + 'center.prof_ave','w')
cout.write('#Average from t=' + str(t[tind])+ ' to t=' + str(t[nt-1]) + '\n')
cout.write('#')
cout.write('[1] x       ')
cout.write('[2] pflux   ')
cout.write('[3] vflux   ')
cout.write('[4] qflux   ')
cout.write('\n')

for i in range (0, nakxc):
  cout.write('%e ' % (dxc*i-0.5*Lxc))
  cout.write('%e ' % pcave[i])
  cout.write('%e ' % vcave[i])
  cout.write('%e ' % qcave[i])
  cout.write('\n')

cout.close()
cout = open(basedir + 'left.prof_ave','w')
cout.write('#Average from t=' + str(t[tind])+ ' to t=' + str(t[nt-1]) + '\n')
cout.write('#')
cout.write('[1] x       ')
cout.write('[2] pflux   ')
cout.write('[3] vflux   ')
cout.write('[4] qflux   ')
cout.write('\n')

for i in range (0, nakxl):
  cout.write('%e ' % (dxc*i-0.5*Lxc))
  cout.write('%e ' % plave[i])
  cout.write('%e ' % vlave[i])
  cout.write('%e ' % qlave[i])
  cout.write('\n')

cout.close()
cout.close()
cout = open(basedir + 'right.prof_ave','w')
cout.write('#Average from t=' + str(t[tind])+ ' to t=' + str(t[nt-1]) + '\n')
cout.write('#')
cout.write('[1] x       ')
cout.write('[2] pflux   ')
cout.write('[3] vflux   ')
cout.write('[4] qflux   ')
cout.write('\n')

for i in range (0, nakxr):
  cout.write('%e ' % (dxc*i-0.5*Lxc))
  cout.write('%e ' % prave[i])
  cout.write('%e ' % vrave[i])
  cout.write('%e ' % qrave[i])
  cout.write('\n')

cout.close()
exit()