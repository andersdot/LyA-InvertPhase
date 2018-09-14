import numpy as np
import sys,os


root = '/home/landerson/src/LyA-InvertPhase/20Mpc_256/'
################################## INPUT #######################################
realizations = 50

redshifts = [2,3,4]

prefix = 'Pk_cdm' #'MF','Pk'
################################################################################

fout_folder = root+'/mean/'
if not(os.path.exists(fout_folder)):  os.system('mkdir '+fout_folder)


for z in redshifts:

    Pk,Pk_NCV = [], []
    for i in xrange(realizations):
        k, Pk_real = np.loadtxt(root+'%d/%s_z=%d.txt'%(i,prefix,z), 
                                unpack=True)
        Pk.append(Pk_real)

        k, Pk_real1 = np.loadtxt(root+'NCV_0_%d/%s_z=%d.txt'%(i,prefix,z), 
                                 unpack=True)

        k, Pk_real2 = np.loadtxt(root+'NCV_1_%d/%s_z=%d.txt'%(i,prefix,z), 
                                 unpack=True)
        Pk_NCV.append(0.5*(Pk_real1 + Pk_real2))

    Pk = np.array(Pk);  Pk_NCV = np.array(Pk_NCV)

    f = open(fout_folder+'%s_mean_z=%d.txt'%(prefix,z), 'w')
    g = open(fout_folder+'%s_mean_NCV_z=%d.txt'%(prefix,z), 'w')
    for j in xrange(k.shape[0]):
        f.write(str(k[j])+' '+str(np.mean(Pk[:,j]))+' '+str(np.std(Pk[:,j]))+\
                '\n')
        g.write(str(k[j])+' '+str(np.mean(Pk_NCV[:,j]))+' '+str(np.std(Pk_NCV[:,j]))+\
                '\n')
    f.close();  g.close()




root = '/home/landerson/src/LyA-InvertPhase/40Mpc_512/'
################################## INPUT #######################################                                                                    
realizations = 50

redshifts = [2,3,4]

prefix = 'Pk_cdm' #'MF','Pk'                                                                                                                        
################################################################################                                                                    

fout_folder = root+'/mean/'
if not(os.path.exists(fout_folder)):  os.system('mkdir '+fout_folder)


for z in redshifts:

    Pk,Pk_NCV = [], []
    for i in xrange(realizations):
        k, Pk_real = np.loadtxt(root+'%d/%s_z=%d.txt'%(i,prefix,z),
                                unpack=True)
        Pk.append(Pk_real)
        if i < 25:
            k, Pk_real1 = np.loadtxt(root+'NCV_0_%d/%s_z=%d.txt'%(i,prefix,z),
                                 unpack=True)

            k, Pk_real2 = np.loadtxt(root+'NCV_1_%d/%s_z=%d.txt'%(i,prefix,z),
                                 unpack=True)
            Pk_NCV.append(0.5*(Pk_real1 + Pk_real2))

    Pk = np.array(Pk);  Pk_NCV = np.array(Pk_NCV)

    f = open(fout_folder+'%s_mean_z=%d.txt'%(prefix,z), 'w')
    g = open(fout_folder+'%s_mean_NCV_z=%d.txt'%(prefix,z), 'w')
    for j in xrange(k.shape[0]):
        f.write(str(k[j])+' '+str(np.mean(Pk[:,j]))+' '+str(np.std(Pk[:,j]))+\
                '\n')
        g.write(str(k[j])+' '+str(np.mean(Pk_NCV[:,j]))+' '+str(np.std(Pk_NCV[:,j]))+\
                '\n')
    f.close();  g.close()
