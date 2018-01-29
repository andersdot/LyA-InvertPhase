print('starting script')

import sys
import astropy.units as u

snapshot_dir = sys.argv[1]
grid_width = 200
spectral_res = 10*u.km/u.s

snapshot_dir_pre = '/mnt/cephtest/landerson/'
print('starting loop')

import PowerSpectraPipeline as psp


for snapshot_num in [4, 7, 9, 12, 14]:
  print(snapshot_dir_pre+snapshot_dir, snapshot_num, grid_width, spectral_res)

  psp.get1dps(snapshot_dir = snapshot_dir_pre + snapshot_dir, snapshot_num=snapshot_num, grid_width=grid_width, spectral_res=spectral_res, reload_snapshot=True, label=None)

  psp.get3dps(snapshot_dir_pre + snapshot_dir, snapshot_num)
