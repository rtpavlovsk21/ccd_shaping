import ccd_utils;
import numpy as np;
import h5py;

def get_ccd_data(file):
    '''
    File is just h5 file with one database. Limits are the adu limits to
    constrain histogram.
    '''
    data = h5py.File(file);
    if( len(data) > 1 ):
        print "warning, only reporting first dbase";
    data = data.items()[0][1];
    print data;
    return data;

def histogram_data(data,limits=[-2**16,2**16]):
    mask = (data > limits[0]) & (data < limits[1]);
    bins = int((data[mask].max()-data[mask].min())/2.);
    n,x = np.histogram(data[mask], bins=bins);
    return n,x;

def selector(file,shift=0):
    '''
    Readin the index vs temperature, return indicies.
    '''
    data = np.genfromtxt(file);
    print data;
    temp = (data[:-1,1]+data[1:,1])/2.;
    count,dtemp=data[:-1,0],np.diff(data[:,1]);
    mask = np.abs(dtemp)<0.01; #0.01Kelvin
    count-=shift;
    mask *= count>0;
    uniq,inds=np.unique(count[mask],return_index=True);
    count,temp = count[mask][inds],np.round(temp[mask][inds]);
    ucount =[];
    utemp = np.unique(temp);
    for t in utemp:
        k = t==temp;
        ucount.append(count[k]);
    return ucount,utemp;

def main(args):
    import matplotlib.pyplot as plt;
    data = get_ccd_data(file=args.file);
    imagenum,temperature=selector(file=args.tempfile)
    for els in imagenum:
        print els;
        n,x = histogram_data(data=data[:,:,els]);
        fit=ccd_utils.gfit_histogram(n,x,0,50);
        print fit;
        plt.semilogy(x[:-1],n);
    plt.show();

if __name__=="__main__":
    import argparse;
    import sys;
    desc = 'Some tools to make plots with ccd data';
    parser = argparse.ArgumentParser(description=desc);
    parser.add_argument('--file',dest='file',default=sys.argv[1],
                        help='H5py file to be processed');
    parser.add_argument('--tempfile',dest='tempfile',default=sys.argv[1],
                        help='Temperature data. Image num vs temp data.');
    args = parser.parse_args();
    main(args=args)
