import sys
import os
import pandas as pd
import numpy as np
import matplotlib as plt
import seaborn as sns
import easygui as eg
import warnings
from difflib import SequenceMatcher
from scipy.stats import ttest_ind

# To compare animals names to see if multiple animals in DB is due to typo
def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

# Turn off annoying known warnings
# Future Warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
# Pandas assign to copy warning
pd.options.mode.chained_assignment = None

# To forward pritn output to log file
class Logger(object):
    def __init__(self,log_file):
        self.terminal = sys.stdout
        self.log = open(log_file, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass    

# To get PPI %
def calcPPI(row,df):
   # print(row)
    if row['Trial_Type']=='PPI':
        stimDB = row['Stim_dB']
        stimMean = df.query('Trial_Type=="Startle" & Stim_dB==@stimDB')['mV Max'].mean()
        return 100*(1-row['mV Max']/stimMean)
    else:
        return row['% PPI']
    
# To get PPI Metrics table
def get_ppi_metrics(x):
    d = {}
    d['Mean Startle (mV)'] = x['mV Max'].mean()
    d['Startle SD (mV)'] = x['mV Max'].std()
    d['Mean PPI (%)'] = x['% PPI'].mean()
    d['PPI SD (%)'] = x['% PPI'].std()
    return pd.Series(d,index=['Mean Startle (mV)',
                              'Startle SD (mV)','Mean PPI (%)',
                              'PPI SD (%)'])

# Set custom plot style
custom_style = {'figure.facecolor':'.8',
                "axes.facecolor":".8",
                'axes.edgecolor':'.8',
                "axes.labelcolor":"black",
                "axes.grid":True,
                'grid.color':'black',
                "text.color":"black",
                "patch_edgecolor":'black',
                "xtick.color":"black",
                "ytick.color":"black",
                'axes.edgecolor':'black'}
sns.set_style('whitegrid',rc=custom_style)
sns.set_context('talk')
col_palette = sns.color_palette("bright")

# Check if command line file input or open file chooser box
if len(sys.argv) < 2:
    file_paths = eg.fileopenbox(msg='Choose PPI datafile to open',filetypes=['*.txt'],multiple=True)
else:
    file_paths = sys.argv[1:]


header= {0:'Chamber',1:'Subject',2:'Session',
         3:'Channel',4:'Trial',5:'Trial Num',
         6:'Group',7:'Param',8:'Trial List Block',
         9:'Samples',10:'Rate',11:'V Start',12:'mV Max',
         13:'T Max',14:'mv Avg',15:'V Peak',16:'T Peak',
         17:'Run Time',18:'TimeStamp',19:'Run Data'}

for fn in file_paths:
    data = pd.read_csv(fn,header=None)
    data = data.rename(columns=header)
    data = data.drop(columns=['Param','Run Data'])
    Nnan = data.isnull().any(axis=1).sum()
    data = data.dropna()
    data_dir = os.path.split(fn)[0]
    
    # Parse Session and Trial 
    # into Date, Animal, Session_Type, Trial_Type, Stim DB, Prepulse DB
    data = pd.concat((data,data.Session.str.split('_', expand=True)),
                     axis=1).rename(columns={0:'Date',1:'Animal',
                                             2:'Session_Type'})
    data = pd.concat((data,data.Trial.str.split('_',expand=True)),
                     axis=1).rename(columns={0:'Trial_Type',1:'Stim_dB',
                                             2:'Prepulse_dB'})
    cols = ['Stim_dB','Prepulse_dB']
    data[cols] = data[cols].apply(pd.to_numeric)
    
    # Print info about data
    animal = data['Subject'].unique()
    if len(animal)>1:
        similarity = [similar(x,y) for i,x in enumerate(animal) for j,y in enumerate(animal) if j>i]
        if all(x>0.9 for x in similarity):
            log_file = os.path.join(data_dir,animal[0]+'_Analysis_Log.txt')
            original_output = sys.stdout
            sys.stdout = Logger(log_file)
            animList = ''.join(x+',' for x in animal)
            print('Found multiple animals names: '+animList)
            print('Similarity was found to be >90% for all names')
            print('Assuming same animal. Changing all names to '+animal[0])
            animal = animal[0]
            data['Subject'] = animal
        else:
            sys.exit('Multiple Animals found in database. Names too dissimilar to assume typo. Quitting.') 
    else:
        log_file = os.path.join(data_dir,animal[0]+'_Analysis_Log.txt')
        original_output = sys.stdout
        sys.stdout = Logger(log_file)
        animal = animal[0]
        
    trial_counts = data.groupby(['Session_Type',
                                'Trial_Type','Stim_dB',
                                 'Prepulse_dB']).size()
    print('Loaded data for '+animal)
    if Nnan >0:
        print(str(Nnan)+' rows were found with unexpected NaN values. These rows were dropped from analysis.')
    print('')
    print('')
    print('Trial Counts')
    print('---------------')
    trial_counts = trial_counts.to_frame().rename(columns={0:'Counts'})
    print(trial_counts.to_string())
    
    # Grab ASR and PPI Data
    asrData = data.query('Session_Type=="ASR"')
    ppiData = data.query('Session_Type=="rPPI"')
    
    if ppiData.shape[0]==0:
        print('')
        print('No rPPI Data found for '+animal)
        noPPI = True
    else:
        noPPI = False
        
        # Calculate PPI 
        ppiData.loc[:,'% PPI'] = np.nan
        ppiData.loc[:,'% PPI'] = ppiData.apply(lambda x: calcPPI(x,ppiData),axis=1)
    
        # Make and print table 
    
        ppi_metrics = ppiData.groupby(['Trial_Type',
                                       'Stim_dB',
                                       'Prepulse_dB']).apply(get_ppi_metrics).round(2)
        print('')
        print('')
        print('PPI Metrics (from rPPI Session)')
        print('------------')
        print(ppi_metrics.to_string())
        
        # Get order for bar plot
        a = ppiData['Trial'].unique()
        sortDF = pd.DataFrame(a,columns=['Trial'])
        sortDF = pd.concat((sortDF,sortDF.Trial.str.split('_',expand=True)),axis=1).rename(columns={0:'Trial_Type',1:'Stim_dB',2:'Prepulse_dB'})
        cols = ['Stim_dB','Prepulse_dB']
        sortDF[cols] = sortDF[cols].apply(pd.to_numeric)
        sortDF['Trial_Type'] = pd.Categorical(sortDF['Trial_Type'],['NoStim','PPIO','Startle','PPI'])
        sortDF = sortDF.sort_values(by=['Trial_Type','Stim_dB','Prepulse_dB'])
        sortDF.reset_index(inplace=True,drop=True)
        sortDF.reset_index(inplace=True)
        sortDF.rename(columns={'index':'x'},inplace=True)
        barOrder = sortDF['Trial']
        
        # Run Statistics
        # Stats: Welch's T-test between PPI Trials and Startle Trials
        ppiStim = ppiData.query("Trial_Type=='PPI'")['Stim_dB'].unique()
        stimDB = [x for x in ppiStim if not ppiData.query('Trial_Type=="Startle" and Stim_dB==@x').empty]
        ppDB = ppiData.query('Trial_Type=="PPI"')['Prepulse_dB'].unique()
        ppiStats = pd.DataFrame(columns=['Stim_dB','Prepulse_dB','T-statistic','p-Value','sigstars','x'])
        for stim in stimDB:
            for pp in ppDB:
                a = ppiData.query('Trial_Type=="Startle" and Stim_dB==@stim')['mV Max']
                b = ppiData.query('Trial_Type=="PPI" and Stim_dB==@stim and Prepulse_dB==@pp')['mV Max']
                if a.empty or b.empty:
                    continue
                res = ttest_ind(a,b,equal_var=False)
                ss = ''
                if res.pvalue<=0.05:
                    ss = '*'
                if res.pvalue<=0.01:
                    ss = '**'
                if res.pvalue<=0.001:
                    ss='***'
                x = sortDF.query('Stim_dB==@stim and Prepulse_dB==@pp')['x'].tolist()[0]
                ppiStats.loc[-1] = [stim,pp,res.statistic,res.pvalue,ss,x]
                ppiStats.reset_index(drop=True,inplace=True)
        
        print('')
        print('')
        print('PPi Statistics (from rPPI Session)')
        print('------------')
        print(ppiStats.drop(columns=['sigstars','x']).to_string())
    
    # Make Figures
    pal = sns.color_palette('bright')
    plt.rcParams["figure.figsize"] = (30, 20)
    plt.rcParams["xtick.labelsize"] = 12
    
    # ASR Plot
    ax = plt.pyplot.subplot(2,1,1)
    g = sns.pointplot(ax=ax,x="Stim_dB",y='mV Max',
                      data=asrData,color='mediumvioletred',ci='sd')
    g.set_title(animal+' Acoustic Startle Response')
    g.set(xlabel='Stimulus dB',ylabel='Max Startle (mV)')

    if noPPI:
        ax = plt.pyplot.subplot(2,1,2)
        ax.axis('off')
        ax.text(0.5,0.5,'No rPPI Data found for '+animal,size=24,ha='center',va='center')
    else:
        # Raw PPI Startle Amplitudes
        ax = plt.pyplot.subplot(2,2,3)
        g = sns.barplot(ax=ax,x='Trial',y='mV Max',data=ppiData,palette='bright',order=barOrder)
        g.set_title(animal+' rPPI Startle Responses')
        g.set(xlabel='Stimulus dB',ylabel='Max Startle (mV)')
        ylim = g.get_ylim()
        errorbars = g.get_lines()
        barTops = [x.get_ydata()[1] for x in errorbars]
        startX = sortDF.x[sortDF.Trial_Type=='Startle'].tolist()[0]
        endX = ppiStats.x[ppiStats['p-Value']<=0.05].tolist()
        stars = ppiStats.sigstars[ppiStats['p-Value']<=0.05].tolist()
        startY = barTops[startX]+10
        endY = [barTops[x] for x in endX]
        for x,y,ss in zip(endX,endY,stars):
            midpoint = (x+startX)/2
            midY = (startY+y)/2
            plt.pyplot.plot([startX,startX,x,x],[startY,startY+5,startY+5,midY],linewidth=2,color='k')
            plt.pyplot.text(midpoint,startY+3,ss)
            startY = startY+12
    
        # % PPI plot
        ax=plt.pyplot.subplot(2,2,4)
        g=sns.barplot(ax = ax,x='Prepulse_dB',y="% PPI",
                      data=ppiData,order = [70,75,80],palette=pal[3:6])
        g.set_title(animal+' Percent PPI')
        g.set(xlabel='Prepulse dB',ylabel='% PPI')

    savefile = os.path.join(data_dir,animal+'_PPI_Results_Figure.png')
    plt.pyplot.savefig(savefile)
    plt.pyplot.close()
    
    print('')
    print('Plot saved to '+savefile)
    print('Done!')
    sys.stdout = original_output
