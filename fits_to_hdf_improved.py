import pyfits;
import numpy as np;
import h5py;
import os;
import time;
import multiprocessing;
from joblib import Parallel, delayed;
import matplotlib.pyplot as plt;

def get_image_array(files):
    header= pyfits.getheader(files[0]);
    numpx = np.asarray([header['naxis2'],header['naxis1']]);
    num_divs=4;
    
    t0=time.time();
    image_array=np.zeros((numpx[0]/2,numpx[1]/2,len(files)),np.int32);
    k=0;
    for fil in files:
        image_array[:,:,k]=pyfits.getdata(fil)[numpx[0]/2:,numpx[1]/2:];
        k+=1;
    return image_array;

def gaus1d(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2));
    
def gfit(vals, lims=[0,0]):
    x=[];n=[];
    mean=[];sigma=[];
    l=[];
    if(lims[0]==lims[1]):
        x,n=np.histogram( vals, bins=( (vals).max()-(vals).min() )/2 );
        maxx=np.argmax(x);
        mean= (n[maxx-1]+n[maxx])/2;
        sigma=50;
        l=np.logical_and( np.less(mean-3*sigma,vals), np.greater(mean+3*sigma,vals) );
    else:
        l=np.logical_and( np.less(lims[0],vals), np.greater(lims[1],vals) );
        mean=np.mean(vals[l]);
        sigma=np.std(vals[l]);
        x,n=np.histogram( vals[l], bins=( (vals[l]).max()-(vals[l]).min() )/2 );
        
    centers=np.asarray([ (n[k-1]+n[k])/2. for k in range(1,len(n)) ]);
    from scipy.optimize import curve_fit
    from scipy import asarray as ar,exp
    popt,pcov = curve_fit(gaus1d,centers,x,p0=[1,mean,sigma])
    return popt;

def readin_files(fits_dir,pre_proc_name):
    files=[];
    if(os.path.exists( os.path.join(fits_dir, pre_proc_name) ) ):
        files=np.genfromtxt( os.path.join(fits_dir,pre_proc_name), dtype='<S50' ) ;
    else:
        files = os.listdir(fits_dir);
        files = filter(os.path.isfile, [ os.path.join(fits_dir,fil) for fil in files ] );
        files = np.asarray(sorted( files,key=lambda x :int(x.split('_')[-1].split('.')[0] ) ));
        fil_lst=file(os.path.join(fits_dir,pre_proc_name) ,'w+');
        print files.dtype;
        for fil in files:
            fil_lst.write(fil+'\n');
        print files[-1];
    return files;

def process(vals):
    return round(gfit(vals)[1]);

fits_dir='raw';
pre_proc_name='file_lst.dat';
files=readin_files(fits_dir,pre_proc_name);

header= pyfits.getheader(files[0]);
numpx = np.asarray([header['naxis2'],header['naxis1']]);

num_cores = multiprocessing.cpu_count()-2;

step=59;
hdf_fil=h5py.File('ccd_images_med_sub.h5','w');
dset=hdf_fil.create_dataset("ba133",(numpx[0]/2,numpx[1]/2,len(files)),dtype=np.int32,maxshape=(numpx[0]/2,numpx[1]/2,None),chunks=(numpx[0]/2,numpx[1]/2,step),compression="lzf")
med=np.median(get_image_array(files[0:step]),axis=2);
sum_=0;
for k in range(step,len(files),step):
    image_array=get_image_array(files[k:k+step]);
    sum_+=(image_array[20:-20,20:-20,:]>2**16-10).sum();
    image_array_prime=image_array-med[:,:,None];

    m=Parallel(n_jobs=num_cores)(delayed(process)(image_array_prime[:,:,i]) for i in range(0,step));
    m=np.asarray(m);

    image_array_prime[ image_array>2**16-10 ]=-1E7;
    print 'Num of pixels at max %s'%sum_;
    dset[:,:,k-step:k]=image_array_prime-m;
    med=np.median(image_array,axis=2);
    print k+step;
    
    
