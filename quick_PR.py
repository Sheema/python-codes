#!/usr/bin/python

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import numpy as np
from scipy import *
from pylab import *
from numpy import *
from scipy import stats
from scipy.optimize import leastsq
from photon_tools import timetag_parse as ttp
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reading the data

donor=ttp.get_strobe_events('/users/sheema/Documents/School/research/data/FRET/2012-12-19/2012-12-19-run_002.timetag',0x1)['t']
acceptor=ttp.get_strobe_events('/users/sheema/Documents/School/research/data/FRET/2012-12-19/2012-12-19-run_002.timetag',0x2)['t']

D_time_events=donor[1024:]/(128.0*1e6) #time in seconds
A_time_events=acceptor[1024:]/(128.0*1e6) #time in seconds

min_D=min(D_time_events)
max_D=max(D_time_events)
min_A=min(A_time_events)
max_A=max(A_time_events)
startAD=max(min_D,min_A)
stopAD=min(max_D,max_A)

""" Histogrammin the arrival times in 5msec bins"""
def raw_data(times,start,stop):
	bin_size=arange(start,stop,0.005)
	counts,bin_edges=histogram(times,bins=bin_size) #5msec bins
	x=bin_edges[:-1]+diff(bin_edges)[0]/2
	return counts, x
	
""" Counts_D are the number of events in each 5msec with a midpoint of xD in Donor channel similar for
     Counts_A and xA for acceptor"""
counts_D,xD=raw_data(D_time_events,startAD,stopAD)
counts_A,xA=raw_data(A_time_events,startAD,stopAD)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#finding average for Background and substracting BG from actual counts

donor_BG=ttp.get_strobe_events('/users/sheema/Documents/School/research/data/FRET/2012-12-19/2012-12-19-run_005.timetag',0x1)['t']
acceptor_BG=ttp.get_strobe_events('/users/sheema/Documents/School/research/data/FRET/2012-12-19/2012-12-19-run_005.timetag',0x2)['t']

D_time_events_BG=donor_BG[1024:]/(128.0*1e6) #time in seconds
A_time_events_BG=acceptor_BG[1024:]/(128.0*1e6) #time in seconds

minBG_D=min(D_time_events_BG)
maxBG_D=max(D_time_events_BG)
minBG_A=min(A_time_events_BG)
maxBG_A=max(A_time_events_BG)
startBG=min(minBG_D,minBG_A)
stopBG=max(maxBG_D,maxBG_A)

BG_D,tD=raw_data(D_time_events_BG,startBG,stopBG)
BG_A,tA=raw_data(A_time_events_BG,startBG,stopBG)

""" mean_BG_D mean value for BG counts in Donor channel similar for mean_BG_A in acceptor channel"""

mean_BG_D=BG_D.mean()
mean_BG_A=BG_A.mean()

#countsnoBG_D and countsnoBG_A are photon counts in Donor and acceptor channel (in 5mese bins) with background subtracted from it.
countsnoBG_D=counts_D-mean_BG_D
countsnoBG_A=counts_A-mean_BG_A


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#setting the threshold and finding the bursts

sum_mean=mean_BG_D+mean_BG_A
counts_sum=countsnoBG_D+countsnoBG_A

threshold1=sum_mean*7
threshold11=sum_mean*8.1
threshold2=sum_mean*9
threshold3=sum_mean*10
threshold4=sum_mean*13
th=array([threshold1,threshold11,threshold2,threshold3,threshold4])

j=0
m=1

counts_list=zeros(len(counts_sum))
index_list=zeros(len(counts_sum))

for i in range(len(counts_sum)):
	if counts_sum[i] >= threshold11:
		counts_list[j]=counts_sum[i]
		index_list[j]=i
		j=j+1
non_zero_counts=counts_list[np.nonzero(counts_list)]
non_zero_index=index_list[np.nonzero(index_list)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#finding the proximity ratio
n=0
PR=zeros(len(non_zero_index))
for k in non_zero_index:
	PR[n]=countsnoBG_A[k]/counts_sum[k]
	n=n+1

min_PR=min(PR)
max_PR=max(PR)
print "Threshold is:", threshold11

bin_size_PR=arange(min_PR,max_PR,0.05)

#savetxt('/Users/sheema/Documents/School/research/Data/FRET/2013-01-03/txt_files/130103_005_T18.txt',PR)

#Plot proximity ratio histograms
PR_counts,PR_values=histogram(PR,bins=bin_size_PR)
hist(PR,bins=bin_size_PR)
xpr=PR_values[:-1]+diff(PR_values)[0]/2
binsize=diff(PR_values)[0]
Area=sum(binsize*PR_counts)
print 'Area=', Area
title('R1R2 200pM- Donly, T=18, 2012-12-19-000')
xlabel('Proximity Ratio',fontsize=19)
ylabel('Number of events',fontsize=19)
xlim(-0.1,1)
ylim(-0.7,0)
legend()
show()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~










