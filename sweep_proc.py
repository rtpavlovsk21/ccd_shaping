from ROOT import TH1D;
import pyfits;
import numpy as np;
import os;
import time;
import cv2;
from matplotlib.pyplot import *;
import multiprocessing;
import time;

sz_lst=[1];
kern_name=['flat','ellipse','cross'];
#kern_name=['flat'];

def make_kernel_list():
    ret_lst=[];
    #sweep flat kernel size
    sz_lst=[1];
    for k in sz_lst:
        ret_lst.append(np.ones((k,k),dtype=np.uint8));
    for k in sz_lst:
        ret_lst.append(cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(k,k)) );
    for k in sz_lst:
        ret_lst.append(cv2.getStructuringElement(cv2.MORPH_CROSS,(k,k)));
    return ret_lst;

def spectra_proc_from_file_list(files,nframe_per_med=29,thread_num=0,semaphore=multiprocessing.Semaphore()):
    header= pyfits.getheader(files[0]);
    numpx = np.asarray([header['naxis2'],header['naxis1']]);
    num_divs=4;
    
    t0=time.time();
    image_array=np.zeros((numpx[0]/2,numpx[1]/2,len(files)),np.int32);
    k=0;
    semaphore.acquire();
    print "(thread %s) Locked hdd access"%thread_num;
    for fil in files:
        image_array[:,:,k]=pyfits.getdata(fil)[numpx[0]/2:,numpx[1]/2:];
        k+=1;
    print "(thread %s) Frame Rate on Readin: %s"%(thread_num,float(len(files))/(time.time()-t0));
    semaphore.release();
    print "(thread %s) Unlocked hdd access"%thread_num;

    nframe_per_med=29;
    nframes=image_array.shape[2]-image_array.shape[2]%nframe_per_med;
    for k in range(0,nframes,nframe_per_med):
        med=np.median(image_array[:,:,k:k+nframe_per_med],axis=2);
        image_array[:,:,k:k+nframe_per_med]-=med[:,:,None]; 
    print "(thread %s) Frame Rate on Med: %s, %s"%(thread_num,float(nframes)/(time.time()-t0),float(nframes));
    
    #cct.make_histogram(image_array=image_array,dims=numpx);
    
    bw_image=(image_array>80).astype(np.uint8);

    kernel_lst=make_kernel_list();
    hist_lst=[];

    kcounter=0;
    for kernel in kernel_lst:
        evt_lst=[];
        nevents=0;
        #kernel=np.ones((3,3),np.uint8);
        for frame in range(0,nframes):
            tbw_image=bw_image[:,:,frame];
            kernel_erode=np.ones((3,3),dtype=np.uint8);
            tbw_image=cv2.erode(tbw_image,kernel_erode,iterations=1);
            tbw_image=cv2.dilate(tbw_image,kernel,iterations=2);
        
            nmarkers,markers=cv2.connectedComponents(tbw_image);
            #marker 0 is probably background
            events=np.zeros(nmarkers-1);
            nevents+=nmarkers-1;
            for k in range(1,nmarkers):
                events[k-1]=np.sum( image_array[markers==k,frame] );
            evt_lst.append(events);
        
        events=np.zeros(nevents);
        for evt in evt_lst:
            events[nevents-evt.shape[0]:nevents]=evt;
            nevents-=evt.shape[0];

        minn=int(0);
        maxx=int(2**20);
        print "(thread %s) Minn: %s Maxx: %s"%(thread_num,minn,maxx);
        name_str=str(thread_num)+str(kcounter)+str(0*np.random.rand(1));
        hist=TH1D(name_str,name_str,maxx-minn,minn,maxx);
        hist.FillN(events.shape[0],events,np.ones(events.shape[0]));
        hist_lst.append(hist);
        kcounter+=1;
    return hist_lst;

def threaded_proc(files,thread_num,queque,semaphore):
    print "Number of files in thread %s is %s"%(thread_num,len(files));
    nframe_per_med=29;
    n_frame=nframe_per_med*20;
    step=int(len(files)/n_frame);
    hist_lst=[];
    for i in range(0,len(files),n_frame):
        if(i==0):
            hist_lst=spectra_proc_from_file_list(files[i:i+n_frame],nframe_per_med=nframe_per_med,thread_num=thread_num,semaphore=semaphore);
        else:
            thist_lst=spectra_proc_from_file_list(files[i:i+n_frame],nframe_per_med=nframe_per_med,thread_num=thread_num,semaphore=semaphore);
            for h in range(0,len(hist_lst)):
                hist_lst[h].Add(thist_lst[h]);
    queque.put(hist_lst);

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

def main():
    #semaphore=threading.Semaphore();
    semaphore=multiprocessing.Semaphore();
    queue=multiprocessing.Queue();
    thread_num=4;
    
    semaphore.acquire();
    fits_dir='raw';
    pre_proc_name='file_lst.dat';
    files=readin_files(fits_dir,pre_proc_name);
    semaphore.release();
    print "done with file readin";
    #files=files[:5000];
    n_t=time.time();
    
    nframes_per_thread=int(len(files)/thread_num);
    for k in range(0,thread_num):
        start=k*nframes_per_thread;
        end=(1+k)*nframes_per_thread;
        t=multiprocessing.Process(target=threaded_proc,args=(files[start:end],k,queue,semaphore));
        t.start();
    
    hist_lst=queue.get();
    for k in range(1,thread_num):
        thist_lst=queue.get();
        for h in range(0,len(hist_lst)):
            hist_lst[h].Add(thist_lst[h]);
    
    drawing=False;
    if(drawing):
        for h in range(0,len( hist_lst)):
            hist_lst[h].Rebin(500);
            hist_lst[h].Scale(100/hist_lst[h].GetEntries());
            if(h==0):
                hist_lst[h].Draw();
                hist_lst[h].SetMaximum(20);
            else:
                hist_lst[h].SetLineColor(h+1);
                hist_lst[h].Draw("SAME");
        raw_input();
    
    for name in kern_name:
        if not os.path.exists(name):
            os.makedirs(name);
    
    for k in range(0,len(hist_lst) ):
        name=kern_name[int(k/(len(sz_lst)))]+'/_'+str(sz_lst[k%len(sz_lst)])+'.root';
        hist_lst[k].SetName( str(sz_lst[k%len(sz_lst)] ));
        hist_lst[k].SaveAs(name);


if __name__=="__main__":
    main();
