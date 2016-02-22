import pyfits;
import numpy as np;
import h5py;
import multiprocessing;
from joblib import Parallel, delayed;
from ccd_utils import *;

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

    quadrant=args.quad;
    step=args.num_im_per_med;
    outfile=args.outfile;
    hdf_fil=h5py.File(outfile,'w');
    dset=hdf_fil.create_dataset(args.dbase_name,
                                (numpx[0]/2,numpx[1]/2,len(files)),
                                dtype=np.int32,
                                maxshape=(numpx[0]/2,numpx[1]/2,None),
                                chunks=(numpx[0]/2,numpx[1]/2,step),
                                compression="lzf")

    med=np.median(get_image_array(files[0:step],quadrant),axis=2);
    sum_=0;
    for k in range(step,len(files),step):
        image_array=get_image_array(files[k:k+step],quadrant);
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
    parser.add_argument('--quad', dest='quad', default=0, type=int,
                        help='Quadrant of image array');
    parser.add_argument('--cpus', dest='num_cpu', default=8, type=int,
                        help='Number of cpus to throw at problem');

    args=parser.parse_args();
    main(args=args);
