#!/usr/bin/python

import datetime
import os
from glob import glob

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42  # to edit text in Illustrator
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd

from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client
from obspy.signal import PPSD


start = UTCDateTime("2018-01-01")
# Leaving UTCDateTime() empty means "now":
end = UTCDateTime("2018-04-25")

ms=start.strftime("%b")
me=(end-1).strftime("%b")
if (ms == me):
   monthstamp=str(ms)
else:
   monthstamp=str(ms)+"-"+str(me)

network = "GT"
station = "BDFB"
location = ""
channel = "BHZ"
dataset = "example"
time_zone = "Brazil/Brasilia"
sitedesc = "in National Park (Brasilia, BR)"

#data_provider = "ODC"
data_provider = "IRIS"
logo = None # 'https://upload.wikimedia.org/wikipedia/commons/thumb/4/44/Logo_SED_2014.png/220px-Logo_SED_2014.png'
#bans = {"2020-03-15 00:00":'Restaurants/Bars/Schools closed', 
#        "2020-03-18 12:00":'Non-essential shops closed'}

bans = {"2020-03-12 06:00":'Universities/Schools closed', 
        "2020-03-20 06:00":'Restaurants/Bars/Non-essential shops closed'}
#bans = None


datelist = pd.date_range(start.datetime, end.datetime, freq="D")
c = Client(data_provider)
print(c)

nslc = "{}.{}.{}.{}".format(network, station, location, channel)
# make sure that wildcard characters are not in nslc
nslc = nslc.replace("*", "").replace("?", "")
for day in datelist:
    datestr = day.strftime("%Y-%m-%d")
    fn = "{}_{}_{}.mseed".format(dataset, datestr, nslc)
    print(fn)
    if day != datelist[-1] and os.path.isfile(fn):
        continue
    else:
        print('AQUI')
#        st = c.get_waveforms(network, station, location, channel,UTCDateTime(day)-1801, UTCDateTime(day)+86400+1801, attach_response=True)
#        print(st)
#        st.write(fn,format="MSEED")
#resp = c.get_stations(UTCDateTime(day), network=network, station=station, location=location, channel=channel, level="response")
resp = c.get_stations(UTCDateTime(day), network='GT', station='BDFB', location='', channel='BHZ', level="response")
print(resp)


for day in datelist:
    datestr = day.strftime("%Y-%m-%d")
    fn_in = "{}_{}_{}.mseed".format(dataset, datestr, nslc)
    if day == datelist[-1] :
        continue
    stall = read(fn_in)
    for mseedid in list(set([tr.id for tr in stall])):
        fn_out = "{}_{}_{}.npz".format(dataset, datestr, mseedid)
        if os.path.isfile(fn_out):
            print("%s done already."%fn_out)
            continue
        st = stall.select(id=mseedid)
        st.attach_response(resp)
        ppsd = PPSD(st[0].stats, metadata=resp,
                    ppsd_length=1800, overlap=0.5,
                    period_smoothing_width_octaves=0.025,
                    period_step_octaves=0.0125,
                    period_limits=(0.008, 50),
                    db_bins=(-200, 20, 0.25))
        ppsd.add(st)
        ppsd.save_npz(fn_out[:-4])
        print(st)
        del st, ppsd
    del stall


ppsds = {}
for day in datelist:
    datestr = day.strftime("%Y-%m-%d")
    fn_pattern = "{}_{}_*.npz".format(dataset, datestr)
    for fn in glob(fn_pattern):
        mseedid = fn.replace(".npz", "").split("_")[-1]
        if mseedid not in ppsds:
            ppsds[mseedid] = PPSD.load_npz(fn)#, allow_pickle=True)
        else:
            ppsds[mseedid].add_npz(fn)#, allow_pickle=True)

yymm=start.strftime("%Y-%m")
[ppsd.plot(yymm+"-ppsd-"+monthstamp+".png",max_percentage=5) for mseedid, ppsd in ppsds.items()]
#[ppsd.plot(yymm+"-ppsd-"+monthstamp+".pdf",max_percentage=5) for mseedid, ppsd in ppsds.items()]
#[ppsd.plot(max_percentage=5) for mseedid, ppsd in ppsds.items()]

#[ppsd.plot_temporal(0.10) for mseedid, ppsd in ppsds.items()]

[ppsd.plot_spectrogram(filename=yymm+"-ppsd-spectrogram-"+monthstamp+".png",clim=(-180,-120)) for mseedid, ppsd in ppsds.items()]
#[ppsd.plot_spectrogram(filename=yymm+"-ppsd-spectrogram-"+monthstamp+".pdf",clim=(-180,-120)) for mseedid, ppsd in ppsds.items()]
#[ppsd.plot_spectrogram(clim=(-180,-120)) for mseedid, ppsd in ppsds.items()]

# Define frequency bands of interest:
freqs = [(0.1,1.0),(1.0,20.0),(4.0,14.0),(4.0,20.0)]


def rms(s, f):
    # Parseval: the RMS in time domain is the sqrt of the integral of the power spectrum
    return np.sqrt(np.trapz(s, f))

displacement_RMS = {}
for mseedid, ppsd in ppsds.items():
    per = ppsd.period_bin_centers
    displacement_RMS[mseedid] = []
    for psd in ppsd.psd_values:
        RMS = {}
        for fmin, fmax in freqs:
            ix = np.where((per>=1.0/fmax) & (per<=1.0/fmin))

            # acceleration power spectrum in Hz
            spec = psd.copy()[ix][::-1]
            f = 1.0/per.copy()[ix][::-1]

            # remove NaNs from the list
            valid = np.where(np.isfinite(spec))[0]
            spec = spec[valid]
            f = f[valid]

            w2f = (2.0 * np.pi * f)

            # The acceleration amplitude spectrum (dB to Power! = divide by 10 and not 20!)
            amp = 10.0**(spec/10.) 

            # velocity spectrum (divide by omega**2)
            vamp = amp / w2f**2

            # displacement spectrum (divide by omega**2)
            damp =  vamp / w2f**2

            RMS["%.1f-%.1f"%(fmin, fmax)] = rms(damp, f)

        displacement_RMS[mseedid].append(RMS)
        index = pd.DatetimeIndex([d.datetime for d in ppsd.times_processed])
    displacement_RMS[mseedid] = pd.DataFrame(displacement_RMS[mseedid], index=index)
    print(mseedid," rms done.")


args = {'band':"4.0-14.0",       # might be None or commented ("4.0-14.0" per default) or any of the tupples in freqs
        'time_zone':time_zone,   # required for clockplots
        'sitedesc':sitedesc,     # might be None or commented
        'logo':logo,             # might be None or commented
        'bans':bans,             # might be None or commented
        'save':'./'              # might be None or commented or a path 
       }

import seismosocialdistancing

seismosocialdistancing.plot(displacement_RMS,type='timeseries', **args)

#seismosocialdistancing.plot(displacement_RMS,type='clockplots',**args)

#seismosocialdistancing.plot(displacement_RMS,type='clockmaps',**args)
