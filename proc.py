#!/usr/bin/python

import obstools
import os
#import numpy as np
#from obspy import UTCDateTime, read, read_inventory, Trace, Stream

file_list = []

for path, subdirs, files in os.walk("./"):
     for name in files:
         if name.startswith("bdfb_202004"):
             #print(os.path.join(path, name))
             file_list.append(os.path.join(path, name))

for f in file_list:
#    obstools.split24to01(f,'displ')
    obstools.prepdata(f)
    print(f)

