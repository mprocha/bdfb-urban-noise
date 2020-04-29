def splitdata(file,time):
    """
    Usage: splitdata(file, time(seconds))
    
    This function read a mseed data and split it 
    in serveral parts with the indicated time
    
    """ 
    import obspy
    from obspy import read, UTCDateTime
 

    st=read(file)
    st.merge(method=1,fill_value=0)

    start=st[0].stats['starttime']

    st1h=st.slice(start,start+time-1)
    st1h.write('teste1.mseed',format='MSEED')



def split24to01(dayfile,option):
    """
    Usage: split24to01(dayfile)
    
    This function read a mseed day-file data and split it 
    in 24 files with 01 hour

    option = displ (displacement)
	   = accel (acceleration)
           = veloc (velocity)
           default (velocity)
    
    """ 
    import obspy
    from obspy import read, read_inventory, UTCDateTime
    from os import mkdir
 

    st=read(dayfile)
    st.merge(method=1,fill_value=0)
    st[0].stats.network='GT'
    st[0].stats.channel='HHZ'
    inv=read_inventory("/SDS/dataless/inventory/GT/GT.BDFB.dataless.xml")

    if option == 'displ':
        dirdata='./1hdispl'
    elif option == 'accel':
        dirdata='./1haccel'
    elif option == 'veloc':
        dirdata='./1hveloc'
    else:
        dirdata='./1hfiles'
        

    sta=st[0].stats['station']
    date=st[0].stats['starttime']
    ref=UTCDateTime(date.date)-1

    x=0
    
    try:
        mkdir(dirdata)
    except:
        pass 

    for i in range(1,25):
        x1=ref+x+1
        x2=ref+x+3600
        x=i*3600
    #    print('%02d' % i)
        print('%s   %s' % (str(x1),str(x2)))
        st1h=st.slice(x1,x2)

        if option == 'displ':
            st1h.remove_response(inventory=inv, output="DISP")
            label='displ'
        elif option == 'accel':
            st1h.remove_response(inventory=inv, output="ACC")
            label='accel'
        elif option == 'veloc':
            st1h.remove_response(inventory=inv, output="VEL")
            label='veloc'
        else:
            label=''
        
        hourfile=str(date.year)+'-'+str("%02d" %date.month)+'-'+str("%02d" %date.day)+'-'+str("%02dh00" %(i-1))+'-'+str("%s"%sta)+'-'+str("%s.mseed"%label)
        st1h.write(dirdata+"/"+hourfile,format='MSEED')

def powerspec(datafile,plot):
    from obspy import read, UTCDateTime
    import numpy as np
    from scipy import fftpack
    from matplotlib import pyplot as plt
    import csv

    st=read(datafile)
    sig=st[0].data
    sig_fft = fftpack.fft(sig)
    power = np.abs(sig_fft)
    time_step=st[0].stats.delta    
    sample_freq = fftpack.fftfreq(sig.size, d=time_step)

    cond=sample_freq>=0

    if plot == True:
       plt.figure(figsize=(18, 5))
       plt.plot(sample_freq[cond], power[cond])
       plt.ylabel('power')
       plt.xlabel('Frequency [Hz]')
       plt.xlim([0, 5])
       plt.show() 
    
    fout=datafile.replace(".mseed",".txt")
    #f=open(fout,"w")
    #f.write("%s %s\n" %(sample_freq[cond], power[cond]))

    s=np.copy(sample_freq[cond])
    p=np.copy(power[cond])
    ls=s.tolist()
    lp=p.tolist()
    z=zip(ls,lp)

    with open(fout, 'w') as f:
        writer = csv.writer(f, delimiter=' ')
        writer.writerows(z)

def read_files(file,spliter,header):
    import numpy as np
    global arrays
    def is_number(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    f=open(file)
    if header == 'y':
        dados=f.readlines()[1:]
    else:
        dados=f.readlines()

    lists=[[] for _ in range(len(dados[0].split(spliter)))]

    for linha in dados:
        dd=linha.split(spliter)
        for i in range(0,len(dd)):
            aa=is_number(dd[i])
            if aa == True:
                lists[i].append(float(dd[i]))
            else:
                lists[i].append(dd[i])
            
    arrays=[]
    for i in range(0,len(dd)):
        arrays.append(np.asarray(lists[i]))
    return(arrays)

def meanspec(a,b,maxfreq,plot):
    #[a,b]=obstools.read_files(file='text.csv',header=None,spliter=' ')
    import numpy as np
    import matplotlib.pyplot as plt
    import csv

    inc=0.01
    cond=a<=(maxfreq+inc)
    ac=a[cond]
    bc=b[cond]

    lista=[];listb=[];x=0;n=0;bsum=0

    for i in range(0,len(ac)):                                        
       acr=float('%.2f'%(np.round(ac[i],decimals=2)))               
       if acr == x:      
          bsum=bsum+bc[i]
          n=n+1
       else:          
          print(float('%.2f'%(np.round(acr-inc,decimals=2))),x)
          bmean=bsum/n      
          lista.append(float('%.20f'%(np.round(acr-inc,decimals=2))))
          listb.append(bmean)
          n=1    
          x=acr     
          bsum=bc[i]

    if plot == True:
       plt.figure(figsize=(18, 5))
       plt.plot(lista,listb)
       plt.ylabel('power')
       plt.xlabel('Frequency [Hz]')
       plt.xlim([0, maxfreq])
       plt.show()

    z=zip(lista,listb)
    
    fout='meanspec.txt'

    with open(fout, 'w') as f:
        writer = csv.writer(f, delimiter=' ')
        writer.writerows(z)

    
