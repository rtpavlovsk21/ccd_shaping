Allowing for simple ccd segmentation for different dilation kernels.

If you find this useful please cite the repo or acknowledge me. 

```
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
```

Right now all you have to do is insert your favorite kernel into the kernel list, add a name, and edit both sz_lst's;