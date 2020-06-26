# -*- coding: utf-8 -*-
"""
# ONA for Python - Noise reduction in Black Carbon data | version 02 (June 2020)

#  ONA is the Optimized Noise reduction Averaging (ONA) algorithm developed by Hagler et al. (2011) to reduce the occurence
# negative values of BC data while preserving the significant dynamics trends in the time series. More information about ONA can be found in:

# Hagler GSW, Yelverton TLB, Vedantham R, Hansen ADA, Turner JR (2011) Post-processing Method to Reduce Noise while Preserving High Time Resolution
# in Aethalometer Real-time Black Carbon Data. Aerosol and Air Quality Research 11, 539–546. doi: 10.4209/aaqr.2011.05.0055

# The original algorithm was developed in MatLab code and this script aims to implement it in Python code.

# The script is available at:
# https://github.com/ncanha/ONA-algorithm-in-Python-for-BC-data-noise-reduction

# This program will look for positive (and final) changes in ATN throughout the data. This code
# anticipates that a variable named “BC_ONA” has already been created and has the following
# columns – 1) serial timestamp, 2) original BC concentrations, and 3) original ATN values.

# The example used in this code is regarding BC data obtained by a microAeth® AE51 (from ARTHLABS) 
# in a monitoring campaign in Greeland.
"""
# Upload the analysis packages:
import numpy as np
import matplotlib.pyplot as plot
import csv
import time
import pandas as pd

# Check your working folder by using "pwd" and then "cd 'path'" to enter in the desired folder
# the raw csv datafile should be in the working folder

# Preparing data
# Remark: the input file as the name of 'input_data.csv' (as an example). It can be updated if the user prefers. 
# The input data should have a similar sctructure to the file "2017_BC_ONA_966.csv"
fields = [] 
rows = []

with open('input_data.csv') as f:
    csvreader = csv.reader(f) 
    fields = next(csvreader)
    
    # Extracting each data row one by one 
    for row in csvreader: 
        rows.append(row)
        
data = np.array(rows)
b= []
data_times =[]
for i in range(len(data[:,0])):
    a=pd.to_datetime(data[i,0])
    b.append(a)
    
    # Time as posix
    data_times.append(time.mktime(a.timetuple()))

# First column is datetime:
data_times = np.array(data_times)

# Second column is black carbon data (BC):
bc = data[:,1]

# Third column is attenuation data (atn):
atn = data[:,2]

# Create three column BC_ONA matrix: 
BC_ONA = np.array([data_times,bc,atn])
BC_ONA = BC_ONA.transpose()    

# Create new variable
BC_ONA_ATN=BC_ONA

# Find filter change points
temp = np.zeros([len(BC_ONA)-1,1])
for i in range(len(BC_ONA)-1):
    temp[i,0] = abs(float(BC_ONA_ATN[i+1,2])-float(BC_ONA_ATN[i,2]))
temp1= np.array(np.where(temp[:,0]>30))
if temp1.size != 0 :
    filtchange = np.array([0])
    filtchange = np.append(filtchange,temp1)
    filtchange = np.append(filtchange,len(BC_ONA_ATN))
else:
    filtchange = np.array([[0],[len(BC_ONA_ATN)]])
del i,temp1

# Check the new variable "temp" that contains the ATN changes between BC measures
temp

# Define the ATN default incremental value
delATN=0.05

# Calculate smoothed BC
temp = np.zeros([len(BC_ONA_ATN),4])
temp[:,:-1] = BC_ONA_ATN
BC_ONA_ATN = temp
del temp
BC_ONA_ATN[:,3] = 1
temp = np.array([])
for i in range(len(BC_ONA_ATN[:,2])):
    if BC_ONA_ATN[i,2] != '-':
        temp2 = float(BC_ONA[i,2]) 
        temp = np.append(temp,temp2)  
        del temp2
for k in range(len(filtchange)-1):
    j = filtchange[k] # Set to first point after filter change
    j = int(j)
    for i in range(j,len(BC_ONA_ATN)):
        if j<filtchange[k+1]:
            if i==j:
                des_ind =np.where(temp[j:len(BC_ONA_ATN)]<=(temp[j]+delATN))
                temp2 = des_ind[-1]
                # des_ind = des_ind[-1]
                # if len(des_ind) == 0:
                if len(temp2) != 0:
                    
                    # Calculated smoothed new BC
                    BC_ONA_ATN[j:temp2[-1]+j+1,1]=np.nanmean(BC_ONA_ATN[j:temp2[-1]+j+1,1])      
                    
                    # Calculate averaging period
                    BC_ONA_ATN[j:temp2[-1]+j+1,3]=len(BC_ONA_ATN[j:temp2[-1]+j+1,1])
                    j = j+temp2[-1]+1
                else:
                    j = j+1
del filtchange

# The variable "BC_ONA_ATN" contains now 4 different collumns: "Time" (as posix), "Corrected BC", "ATN" and "points average"

# METRICS OF SMOOTHING PERFORMANCE 
# The code from this point forwards only displays the results and does not do any data alteration.

#1. Reduction of negatives - fraction in original vs. remaining
temp2 = 0
temp3 = 0

for i in range(len(BC_ONA[:,1])):
    if BC_ONA_ATN[i,1] != '':
        if float(BC_ONA[i,1])<0:
            temp2 += 1
    if BC_ONA_ATN[i,1] != '-':
            if float(BC_ONA_ATN[i,1])<0:
             temp3 += 1     

numneg_org=temp2/len(BC_ONA[:,1])
numneg_filt=temp3/len(BC_ONA[:,1])
del temp2,temp3
print("Fraction of negative values in original dataset = ",numneg_org)
print("Fraction of negative values after ONA correction = ",numneg_filt)
del numneg_org,numneg_filt

#2. Reduction of noise
temp = np.zeros([len(BC_ONA)-1,2])
for i in range(len(BC_ONA)-1):
    if BC_ONA_ATN[i,1] != '' and BC_ONA_ATN[i+1,1] != '':
        temp[i,0]=abs(float(BC_ONA[i+1,1])-float(BC_ONA[i,1]))
    if BC_ONA_ATN[i,1] != '-' and BC_ONA_ATN[i+1,1] != '-':
        temp[i,1]=abs(float(BC_ONA_ATN[i+1,1])-float(BC_ONA_ATN[i,1]))

noise = np.array([np.nanmean(temp[:,0]),np.nanmean(temp[:,1])])
langs = ['Original dataset','ONA corrected dataset']
fig = plot.figure()
plot.bar(langs,noise)
plot.title ('Reduction of Noise between BC measurements')
plot.ylabel('Avg Abs BC$_{t+1 - t}$ (ng.$m^-$$^3$)')
plot.rcParams['figure.figsize']=[12,6]
plot.savefig('Figure 1 - Noise Reduction.png')
print('BC Noise reduction (from to, in ng.m-3):', noise)

#3. Averaging interval and histogram of points averaged
time_inc = abs(float(BC_ONA[1,0])-float(BC_ONA[0,0]))
temp = BC_ONA_ATN[:,3].astype(np.float)
timeavg = np.array([time_inc,np.nanmean(temp)*time_inc])
fig2 = plot.figure()
plot.title ('Averaged time after ONA correction')
plot.bar(langs,timeavg)
plot.ylabel('avg timebase (s)')
plot.savefig('Figure 2 - Averaged time after ONA correction.png')

del time_inc,temp,timeavg
fig3 = plot.figure()
plot.hist(BC_ONA_ATN[:,3], bins = 10)
plot.title ('Histogram of points averaged')
plot.xlabel('# pts averaged')
plot.ylabel('N pts affected')
plot.savefig('Figure 3 - Histogram of points averaged.png')

#4. Before - top figure BC time series, bottom figure ATN time series

fig4 = plot.figure()

# posix times and datetimes
plot.subplot(211)
yticks = [-2000,0,2000]
plot.plot(BC_ONA[::10,0], BC_ONA[::10,1],linewidth=0.75)
plot.yticks(yticks)
plot.ylabel('BC (ng.$m^-$$^3$)')
plot.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)

# posix times:
plot.subplot(212)
plot.plot(BC_ONA[::10,0], BC_ONA[::10,2],linewidth=0.75)
plot.ylabel('ATN')
plot.xlabel('Date')
plot.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
plot.savefig('Figure 4 - Before BC and ATN time series.png')

#5 After - top figure BC and BC-ONA time series, bottom figure averaging time
fig5 = plot.figure()
plot.subplot(221)
yticks = [-2000,0,2000]
plot.yticks(yticks)
plot.plot(BC_ONA[::10,0], BC_ONA[::10,1],color ='black',linewidth=0.5)
plot.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=True, labelleft=False)
plot.ylabel('average BC concentration (ng m^-^3)')
plot.axis('equal')
plot.subplot(222)

plot.plot(BC_ONA_ATN[:,0], BC_ONA_ATN[:,1],color ='green', marker='o')
plot.xlabel('time')
plot.axis([min(BC_ONA_ATN[:,0]), max(BC_ONA_ATN[:,0]),0,300])
plot.ylabel('average ONA BC concentration (ng.$m^-$$^3$)')

plot.subplot(223)
plot.plot(BC_ONA_ATN[:,0], BC_ONA_ATN[::,3],color ='blue',linewidth=1)
plot.xlabel('time')
plot.ylabel('# pts averaged')
plot.savefig('Figure 5 - BC and BC-ONA time series.png')

# Conversion of Numpy array to Pandas dataframe of BC_ONA_ATN
Final_BC_ONA = pd.DataFrame({'Timestamp': BC_ONA_ATN[:,0],'BC ONA': BC_ONA_ATN[:,1],'ATN': BC_ONA_ATN[:,2],'Points averaged': BC_ONA_ATN[:,3],})

# Check the number of occurencies of negative concentrations of BC:
neg=Final_BC_ONA[['BC ONA']].agg(lambda x: sum(x < 0)).sum()
print(neg)
# Total number of concentrations values available for BC:
total=len(Final_BC_ONA['BC ONA'])
print(total)
# Calculate the percentage of BC negative values:
pneg = neg/total*100
print("The percentage of negative BC values in the dataset is", pneg.round(0).astype(int), "%." )

# Conversion of timestamp from Posix to Year-Month-Day Hour:Minute
Final_BC_ONA['Timestamp'] = pd.to_datetime(Final_BC_ONA['Timestamp'],unit='s', dayfirst=True )

# Histogram of ONA BC data (range limits of visualisation should be updated taking in account the working dataset)

Hist_BC_ONA=Final_BC_ONA[['BC ONA']]
fig6 = Hist_BC_ONA.plot.hist(range=[-5, 115],bins=[-5, 0,5,10,15, 20,25, 30,35, 40,45, 50,55,60,65,70,75,80,85,90,95, 100, 105, 110, 115], title = 'Histogram of ONA BC data')
fig6.set_xlabel('BC Concentration (ng.$m^-$$^3$)')

# To save the figure
plot.savefig('Figure 6 - Histogram of ONA BC data')


# Check the histogram using a narrow timeframe:
fig7 = Hist_BC_ONA.plot.hist(range=[0, 20], bins=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],title = 'Histogram of selected ONA BC data')
fig7.set_xlabel('BC Concentration (ng.$m^-$$^3$)')

# To save the figure
plot.savefig('Figure 7 - Histogram of selected ONA BC data')

# Goal: compare BC concentrations before and after ONA correction
# 2.1. How to create a new dataframe (BC_ba: before and after) with the timestamp, BC data before and after correction, ATN and number of averaged points:

# 2.1.1. In order to avoid problems with the "timestamp", use the time collumn from the original csv file:
Initial_BC = pd.read_csv('2017_BC_ONA_966.csv', names=['Timestamp','BC','ATN'])

# 2.1.2. Create the new dataframe:
BC_ba = Initial_BC[['Timestamp']].join(Initial_BC[['BC']].join(Final_BC_ONA[['BC ONA']].join(Final_BC_ONA[['ATN']].join(Final_BC_ONA[['Points averaged']]))))

# 2.1.3. "BC" is not in the correct format (float) but it is in the "object" format. To change:
BC_ba["BC"] = BC_ba.BC.astype(float)

# 2.1.4. Check data types:
BC_ba.dtypes

# 3. Graphs of BC before and after of ONA implementation:
fig8 = BC_ba.plot('Timestamp', 'BC', title = 'Raw BC data during June2017 monitoring campaign',
                   ylim=[-2000,2500])
fig8.set_xlabel('Date')
fig8.set_ylabel('BC Concentration (ng.$m^-$$^3$)')
fig8.grid(True)
# To save the figure
plot.savefig('Figure 8 - Raw BC data')

fig9 = BC_ba.plot('Timestamp', 'BC ONA', title = 'Corrected BC data during June2017 monitoring campaign',
                   ylim=[0,120])
fig9.set_xlabel('Date')
fig9.set_ylabel('BC Concentration (ng.$m^-$$^3$)')
fig9.grid(True)
# To save the figure
plot.savefig('Figure 9 - ONA corrected BC data')

# 4. Graph of ATN
fig10 = BC_ba.plot('Timestamp', 'ATN', title = 'ATN during June2017 monitoring campaign',
                   ylim=[4,12.5])
fig10.set_xlabel('Date')
fig10.set_ylabel('Attenuation')
fig10.grid(True)

plot.savefig('Figure 10 - Attenuation during monitoring period')

# 5. Writing and saving file with all data

# This final csv file will have 5 different collumns: Timestamp, BC, BC ONA, ATN and Points averaged:

selection = BC_ba

output_selection = 'BC_ba_ONA_final.csv'

# Save dataframe to csv
selection.to_csv(output_selection, sep=',', index = False)

# "index = False" will remove the column index

# 6. Get descritive statistics summary only for "Temp":
stats = BC_ba[['BC','BC ONA','ATN','Points averaged']].describe()
print(stats)

# At the end of the script, the working folder should have the input file (csv), the output file with all data (csv) and 10 Figures regarding the 
# the performance of the ONA algorithm, raw and corrected data.
# The end!
print('Done! BC data has been corrected using the ONA algorithm for noise reduction!')