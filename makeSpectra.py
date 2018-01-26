import PowerSpectrPipeline as psp
import sys

snapshot_dir = sys.argv[1]
grid_width = 200
spectral_res = 10*u.km/u.s

for snapshot_num in range(15):
  psp.get1dps(snapshot_dir = snapshot_dir, snapshot_num=snapshot_num, grid_width=20, spectral_res=50*u.km/u.s, reload_snapshot=True, label=None)
  psp.get3dps(snapshot_dir, snapshot_num)
