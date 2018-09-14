
import numpy as np
import sys,os
import readgadget
import MAS_library as MASL
import Pk_library as PKL


def generatePk(simroot, outroot, dims, ptypes, MAS, do_RSD, axis, threads, bigbox=True):
    for i in xrange(0, 50):

        for prefix in ["","NCV_0_","NCV_1_"]:
            
            if ((prefix == "NCV_0_") or (prefix == "NCV_1_")) and (i > 24) and bigbox : continue
            
            for num in [0,1,2]:

                print 'Working with snap %s%d'%(prefix,i)

                snapshot = simroot+prefix+'%d/snap_%03d'%(i,num)

                header   = readgadget.header(snapshot)
                BoxSize  = header.boxsize/1e3  #Mpc/h                                                                         
                redshift = header.redshift
                masses   = header.massarr*1e10  #Msun/h                                                          
                fout_dir = outroot+prefix+'%d'%(i)
                fout = fout_dir + '/Pk_cdm_z=%d.txt'%(round(redshift))
                
                if os.path.exists(fout):  continue
                if not os.path.isdir(fout_dir): os.mkdir(fout_dir)
                
                delta_m = MASL.density_field_gadget(snapshot, ptypes, dims, MAS, do_RSD, axis)
                delta_m /= np.mean(delta_m, dtype=np.float64);  delta_m -= 1.0
                
                Pk = PKL.Pk(delta_m, BoxSize, axis, MAS, threads)
                np.savetxt(fout, np.transpose([Pk.k3D, Pk.Pk[:,0]]))
                

################################ INPUT #################################
simroot = '/simons/scratch/fvillaescusa/Lya_ncv/20Mpc_256/'
outroot = '/home/landerson/src/LyA-InvertPhase/20Mpc_256/'

dims    = 256
ptypes  = [1]
MAS     = 'CIC'
do_RSD  = False 
axis    = 0 
threads = 16
########################################################################

generatePk(simroot, outroot, dims, ptypes, MAS, do_RSD, axis, threads, bigbox=False)


################################ INPUT #################################                                                                            
simroot = '/simons/scratch/fvillaescusa/Lya_ncv/40Mpc_512/'
outroot = '/home/landerson/src/LyA-InvertPhase/40Mpc_512/'

dims    = 512
ptypes  = [1]
MAS     = 'CIC'
do_RSD  = False
axis    = 0
threads = 16
########################################################################                                                                            
generatePk(simroot, outroot, dims, ptypes, MAS, do_RSD, axis, threads)

