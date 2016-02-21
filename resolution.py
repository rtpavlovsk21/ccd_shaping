import numpy as np;
from matplotlib.pyplot import *;
import sys;
#args should be energy of peak_ene, isotope

def make_sweep_plots(peak_ene_kev,isotope):
    args=[' ',peak_ene_kev,isotope];
    dir_lst=['cross','flat','ellipse'];
    fig,ax1=subplots();
    data_lst=[];
    for dir in dir_lst:
        data_lst.append(np.genfromtxt(dir+'/resolution.txt',delimiter=','));
        data_lst[-1]=np.sort(data_lst[-1],axis=0);
    print data_lst[0].shape
    k=0;
    for data in data_lst:
        ax1.loglog(data[:,0],data[:,1]*100.,label=dir_lst[k]);
        ax1.set_ylabel('FWHM %');
        ax1.set_xlabel('Kernel Size');
        ax1.set_title('FWHM vs Kernel Size for %skeV Gamma %s'%(args[1],args[2]) )
        k+=1;
    legend()
    savefig('fwhm_vs_kernel_percent_%s_%s.png'%(args[1],args[2].replace(' ','_')),transparent=True);
    figure();
    fig,ax2=subplots();
    k=0;
    for data in data_lst:
    #    ax2.loglog(data[:,0],data[:,1]*356.,label=dir_lst[k]);
        ax2.loglog(data[:,0],data[:,1]*float(args[1]),label=dir_lst[k]);
        ax2.set_ylabel('FWHM (keV)');
        ax2.set_xlabel('Kernel Size');
        ax2.set_title('FWHM vs Kernel Size for %skeV Gamma %s'%(args[1],args[2]) );
        k+=1;
    legend();
    savefig('fwhm_vs_kernel_energy_%s_%s.png'%(args[1],args[2].replace(' ','_')),transparent=True);

def main():
    args=sys.argv;
    make_sweep_plots(args);

if __name__=='__main__':
    main();
