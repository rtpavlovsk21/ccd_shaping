import pyfits;
import numpy as np;
import h5py;
import os;
import time;
import multiprocessing;
from joblib import Parallel, delayed;
import matplotlib.pyplot as plt;
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

def get_image_array(files):
    '''
    Take in files list and return stacked ccd images MxNxlen(files). Note
    that the full ccd image is 2Mx2N, this only extracts one quadrant.
    '''
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
    '''
    Gaussian for fitting
    '''
    return a*np.exp(-(x-x0)**2/(2*sigma**2));

def gfit(vals, lims=[0,0]):
    '''
    Flatten vals, ccd pixel values, and fit them with a guassian for black
    level noise.
    '''
    x=[];n=[];
    mean=[];sigma=[];
    l=[];
    if(lims[0]==lims[1]):
        bins = int( (vals.max()-vals.min())/2. );
        n,x=np.histogram( vals, bins=bins );
        maxn=np.argmax(n);
        mean= (x[maxn-1]+x[maxn])/2.;
        sigma=50;
        l= (mean-3*sigma<vals) & (mean+3*sigma>vals) ;
    else:
        l= (lims[0]<vals) & (lims[1]>vals);
        mean=np.mean(vals[l]);
        sigma=np.std(vals[l]);
        bins = int( (vals[l].max()-vals[l].min())/2. );
        n,x=np.histogram( vals[l], bins=bins);

    centers=np.asarray([ (x[k-1]+x[k])/2. for k in range(1,len(x)) ]);
    popt,pcov = curve_fit(gaus1d,centers,n,p0=[1,mean,sigma])
    return popt;

def readin_files(fits_dir,pre_proc_name):
    '''
    Take in fits_dir to collect all fits with pre_proc_name cached file names
    '''
    files=[];
    if(os.path.exists( os.path.join(fits_dir, pre_proc_name) ) ):
        files=np.genfromtxt(os.path.join(fits_dir,pre_proc_name), dtype='<S60');
    else:
        files = os.listdir(fits_dir);
        files =  [ os.path.join(fits_dir,fil)
                   for fil in files
                   if( fil.endswith('fits') | fil.endswith('fit') )];
        files = filter(os.path.isfile,files);
        files = np.asarray(sorted( files,key=lambda x :int(x.split('_')[-1].split('.')[0] ) ));
        fil_lst=file(os.path.join(fits_dir,pre_proc_name) ,'w+');
        for fil in files:
            fil_lst.write(fil+'\n');
        print files[-1];
    return files;

def process(vals):
    '''
    Fit a gaussian to vals and return mean. This works to subtract dc offset
    after median subtraction.
    '''
    return round(gfit(vals)[1]);

def main(args):

    fits_dir=args.fits_dir;
    pre_proc_name=args.names_list;
    files=readin_files(fits_dir,pre_proc_name);

    header= pyfits.getheader(files[0]);
    numpx = np.asarray([header['naxis2'],header['naxis1']]);

    num_cores = args.num_cpu;

    step=args.num_im_per_med;
    outfile=args.outfile;
    hdf_fil=h5py.File(outfile,'w');
    dset=hdf_fil.create_dataset(args.dbase_name,
                                (numpx[0]/2,numpx[1]/2,len(files)),
                                dtype=np.int32,
                                maxshape=(numpx[0]/2,numpx[1]/2,None),
                                chunks=(numpx[0]/2,numpx[1]/2,step),
                                compression="lzf")

    med=np.median(get_image_array(files[0:step]),axis=2);
    sum_=0;
    for k in range(step,len(files),step):
        image_array=get_image_array(files[k:k+step]);
        sum_+=(image_array[20:-20,20:-20,:]>2**16-10).sum();
        image_array_prime=image_array-med[:,:,None];

        shape = image_array_prime.shape[2];
        m=Parallel(n_jobs=num_cores)(delayed(process)(image_array_prime[:,:,i]) for i in range(0,shape));
        m=np.asarray(m);

        image_array_prime[ image_array>2**16-10 ]=-1E7;
        print 'Num of pixels at max %s'%sum_;
        dset[:,:,k-shape:k]=image_array_prime-m;
        med=np.median(image_array,axis=2);
        print k-step, k;

if __name__ == '__main__':
    import argparse;

    dsc = 'Median subtract fits files, save the result as an hdf5 files'
    parser = argparse.ArgumentParser(description=dsc);
    parser.add_argument('--dir', dest='fits_dir', help='Directory of fits');
    parser.add_argument('--list', dest='names_list',
                        help='File with fits names list');
    parser.add_argument('--dbase', dest='dbase_name', default='ccd_dbase',
                        help='Hdf database name');
    parser.add_argument('--outfile', dest='outfile',
                        default='ccd_sub.h5',help='Output file h5');
    parser.add_argument('--nipm', dest='num_im_per_med', default=59, type=int,
                        help='Number of images per median/chunk');
    parser.add_argument('--cpus', dest='num_cpu', default=8, type=int,
                        help='Number of cpus to throw at problem');

    args=parser.parse_args();
    main(args=args);
