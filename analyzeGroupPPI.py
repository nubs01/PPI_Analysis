import sys, os, itertools, warnings
import pandas as pd
import numpy as np
import matplotlib as plt
import seaborn as sns
import easygui as eg
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
        tmp = df.query('Trial_Type=="Startle" & Stim_dB==@stimDB')
        stimMean = tmp['Norm_Startle'].mean()
        return 100*(1-row['Norm_Startle']/stimMean)
    else:
        return row['% PPI']
    
# To get PPI Metrics table
def get_ppi_metrics(x):
    d = {}
    d['Mean Startle (mV)'] = x['Norm_Startle'].mean()
    d['Startle SD (mV)'] = x['Norm_Startle'].std()
    d['Mean PPI (%)'] = x['% PPI'].mean()
    d['PPI SD (%)'] = x['% PPI'].std()
    return pd.Series(d,index=['Mean Startle (mV)',
                              'Startle SD (mV)','Mean PPI (%)',
                              'PPI SD (%)'])

def get_ppi_data(fn):
    header= {0:'Chamber',1:'Subject',2:'Session',
             3:'Channel',4:'Trial',5:'Trial Num',
             6:'Group',7:'Param',8:'Trial List Block',
             9:'Samples',10:'Rate',11:'V Start',12:'mV Max',
             13:'T Max',14:'mV Avg',15:'V Peak',16:'T Peak',
             17:'Run Time',18:'TimeStamp',19:'Run Data'}
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
    data = pd.concat((data,data.Group.str.split('_',expand=True)),
                    axis=1).rename(columns={0:'Trash',1:'GenGeno'})

    cols = ['Stim_dB','Prepulse_dB']
    data[cols] = data[cols].apply(pd.to_numeric)
    
    # Print info about data
    animal = data['Subject'].unique()
    if len(animal)>1:
        similarity = [similar(x,y) for i,x in enumerate(animal) for j,y in enumerate(animal) if j>i]
        if all(x>0.9 for x in similarity):
            animList = ''.join(x+',' for x in animal)
            print('Found multiple animals names: '+animList)
            print('Similarity was found to be >90% for all names')
            print('Assuming same animal. Changing all names to '+animal[0])
            animal = animal[0]
            data['Subject'] = animal
        else:
            sys.exit('Multiple Animals found in database. Names too dissimilar to assume typo. Quitting.') 
    else:
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
    
    # Normalize Startle response
    data['Norm_Startle'] = data['mV Max']-data['mV Avg']
    
    # Grab ASR and PPI Data
    asrData = data.query('Session_Type=="ASR"')
    ppiData = data.query('Session_Type=="rPPI"')
    if ppiData.shape[0]==0:
        print('')
        print('No rPPI Data found for '+animal)
    else:
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
    return ppiData


def get_ppi_stats(data,groups=None):
    # Run Statistics
    # Stats: Welch's T-test between GenGeno groups for % PPI at each prepulse volume
    ppiDat = data.query('Trial_Type=="PPI"')
    if groups is None:
        groups = ppiDat['GenGeno'].unique().tolist()
    ppDB = [70,75,80] #HARDCODED
    comparisons = list(itertools.combinations(groups,2))
    ppiStats = pd.DataFrame(columns=['Group A','Group B','Prepulse_dB','T-statistic','p-Value','sigstars','xA','xB','ax','dist'])
    for pp in ppDB:
        for test_grps in comparisons:
            a = ppiDat.query('GenGeno==@test_grps[0] and Prepulse_dB==@pp')['% PPI']
            b = ppiDat.query('GenGeno==@test_grps[1] and Prepulse_dB==@pp')['% PPI']
            if a.empty or b.empty:
                continue
            res = ttest_ind(a,b,equal_var=False)
            ss = ''
            if res.pvalue<=0.05:
                ss = '*'
            if res.pvalue<=0.01:
                ss = '**'
            if res.pvalue<=.001:
                ss = '***'
            xA = groups.index(test_grps[0])
            xB = groups.index(test_grps[1])
            ax = ppDB.index(pp)
            dist = abs(xA-xB)
            ppiStats.loc[-1] = [test_grps[0],test_grps[1],pp,res.statistic,res.pvalue,ss,xA,xB,ax,dist]
            ppiStats.reset_index(drop=True,inplace=True)
    ppiStats = ppiStats.sort_values(by=['ax','dist'],ascending=True)
    return ppiStats,groups

if __name__ == '__main__':

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
        file_paths = eg.fileopenbox(msg='Choose all PPI csv files to load',filetypes=['*.txt','*.csv'],multiple=True)
    else:
        file_paths = sys.argv[1:]

    data_dir = os.path.split(file_paths[0])[0]
    log_file = os.path.join(data_dir,'Group_Analysis_Log.txt')
    original_output = sys.stdout
    sys.stdout = Logger(log_file)

    data = None
    for fn in file_paths:
        tmp = get_ppi_data(fn)
        if tmp is not None:
            if data is None:
                data = tmp
            else:
                data = data.append(tmp)

    ppiStats,groups = get_ppi_stats(data)
    print('')
    print('Statistics, T-test')
    print(ppiStats)
    prepulseDB = [70,75,80]



    # Make Figures
    pal = sns.color_palette('bright')
    plt.rcParams["figure.figsize"] = (30, 20)
    plt.rcParams["xtick.labelsize"] = 12
    
    # ASR Plot
    g = sns.catplot(x="GenGeno",y='% PPI',kind='bar',order=groups,
                    data=data,col='Prepulse_dB',col_order=prepulseDB)
    plt.pyplot.subplots_adjust(top=0.8)
    g.fig.suptitle('Percent Prepulse Inhibition')
    axes = g.axes.flatten()
    axes[0].set_title('Prepulse dB: 70')
    axes[1].set_title('Prepulse dB: 75')
    axes[2].set_title('Prepulse dB: 80')
    axes[0].set_ylabel('% PPI')
    axes[0].set_xlabel('')
    axes[1].set_xlabel('Gender x Genotype')
    axes[2].set_xlabel('')

    for i,ax in enumerate(axes,start=0):
        plt.pyplot.sca(ax)
        errorbars = ax.get_lines()
        barTops = [x.get_ydata()[1] for x in errorbars]
        yVals = [x+5 for x in barTops]
        pp = prepulseDB[i]
        ps = ppiStats.query('Prepulse_dB==@pp')
        lineYs = []
        for row in ps.iterrows():
            if row[1]['p-Value'] <=0.05:
                g1 = groups.index(row[1]['Group A'])
                g2 = groups.index(row[1]['Group B'])
                y1 = yVals[g1]
                y2 = yVals[g2]
                yBar = y1+5 if y1>y2 else y2+5
                if any(t<=yBar+5 and t>=yBar-5 for t in lineYs):
                    yBar = max(lineYs)+7
                lineYs.append(yBar)
                yVals[g1]+=7
                yVals[g2]+=7
                midX = (g1+g2)/2
                plt.pyplot.plot([g1,g1,g2,g2],[y1,yBar,yBar,y2],linewidth=2,color='k')
                plt.pyplot.text(midX,yBar+1,row[1]['sigstars'],ha='center',va='center')


    savefile = os.path.join(data_dir,'Group PPI Plot.png')
    plt.pyplot.savefig(savefile)
    plt.pyplot.close()
    print('PPI Plot save to '+savefile)
    
    plt.rcParams["figure.figsize"] = (10, 10)
    g = sns.barplot(data = data.query('Trial_Type=="Startle"'),x='GenGeno',y='Norm_Startle')
    g.set_title('Group Startle Responses')
    g.set(xlabel='Gender & Genotype',ylabel='Startle Amplitude (mV)')
    savefile = os.path.join(data_dir,'Group Startle Response Plot.png')
    plt.pyplot.savefig(savefile)
    plt.pyplot.close()

    print('Startle Response plot saved to '+savefile)
    sys.stdout = original_output

