import pyfits;
import numpy as np;
import h5py;
import matplotlib.pyplot as plt;
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

def quad_limits(quadrant,numpx):
    lims=[];
    if(quadrant==0):
        lims=[0,numpx[0]/2, numpx[1]/2,numpx[1]];
    elif(quadrant==1):
        lims=[0,numpx[0]/2, 0,numpx[1]/2];
    elif(quadrant==2):
        lims=[numpx[0]/2,numpx[0], 0,numpx[1]/2];
    elif(quadrant==3):
        lims=[numpx[0]/2,numpx[0], numpx[1]/2,numpx[1]];
    return np.asarray(lims);

def get_image_array(files,quadrant):
    '''
    Take in files list and return stacked ccd images MxNxlen(files). Note
    that the full ccd image is 2Mx2N, this only extracts one quadrant.

       Quadrant is [0,3] to select quad.
    _____________
    |     |     |
    |  1  |  0  |
    |     |     |
    |-----------|
    |     |     |
    |  2  |  3  |
    |___________|
    '''
    header= pyfits.getheader(files[0]);
    numpx = np.asarray([header['naxis2'],header['naxis1']]);
    num_divs=4;

    q_lims=quad_limits(quadrant,numpx);
    image_array=np.zeros((q_lims[1]-q_lims[0],q_lims[3]-q_lims[2],len(files)),np.int32);

    k=0;
    for fil in files:
        image_array[:,:,k]=pyfits.getdata(fil)[q_lims[0]:q_lims[1],q_lims[2]:q_lims[3]];
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

    return gfit_histogram(n,x,mean,sigma);

def gfit_histogram(n,x,mean_guess,sigma_guess):
    centers=np.asarray([ (x[k-1]+x[k])/2. for k in range(1,len(x)) ]);
    popt,pcov = curve_fit(gaus1d,centers,n,p0=[1,mean_guess,sigma_guess])
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

