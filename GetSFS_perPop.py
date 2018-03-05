##SFS PLOT
import pickle
import matplotlib.pyplot as plt
from operator import itemgetter
import pandas as pd
import glob
from collections import Counter
import sys
import os
import argparse

def main(args):
    names=pd.read_table(
        '/Users/luke/bin/smaller_mut_spectrum_pipeline/Name.Pop.Day.txt', 
        sep=' ', names=['Name','Pop_Day'])
    Pops=list(names['Pop_Day'].unique())

    AF=Counter()
    list_of_files = glob.glob(
        args.inputDir+'Chr*_AlleleFrequencyPerPopDay.txt')

    for file_name in list_of_files:
        if os.path.getsize(file_name) > 0:
            print(file_name)
            if sys.getsizeof(AF) > 800000000:
                print('##Output Dump 1')
                AFname1 = args.outDir +'AlleleFrequencyDUMP1.txt'
                with open(AFname1, 'wb') as outputfile:
                    pickle.dump(AF, outputfile)
                AF = Counter()
            AF = AF + pickle.load(open(file_name, "rb"))
            print(sys.getsizeof(AF))
            
    print('##Output Dump 2')
    AFname2 = args.outDir +'AlleleFrequencyDUMP2.txt'
    with open(AFname2, 'wb') as outputfile:
        pickle.dump(AF, outputfile)
    
    AF = AF + pickle.load(open(args.out +'AlleleFrequencyDUMP1.txt', "rb"))
    D=pd.DataFrame.from_dict(AF,orient='index').reset_index()
    D[Pops] = D['index'].apply(pd.Series)
    D=D.drop(['index'], axis=1)
    D.columns.values[0] = 'Count'

    for i in range(0,52,2):
        print(Pops[i][:3])
        Sub=D.loc[:,['Count',Pops[i], Pops[i+1]]]
        Sub=Sub.groupby([Pops[i], Pops[i+1]]).sum().reset_index()
        Sub.to_csv(args.outDir+'AllChr'
                   +Pops[i][:3]+'HeatMap.txt')
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'GET SFS per POP')
    parser.add_argument('-inputDir', help = 'input dir')
    parser.add_argument('-outDir', help= 'output directory')
    args = parser.parse_args()
    main(args)
