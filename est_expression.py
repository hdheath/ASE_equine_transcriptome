import os
import pandas as pd

def CountUp():
    
    df_list = []  # This will hold DataFrames which we will concatenate at the end
    directory = os.path.expanduser("./ASDEG/results/samtools")
    files = [f for f in os.listdir(directory) if f.endswith('results')]
    
    low_quality_symbols = set([chr(i) for i in range(33, 33+20)])
    
    for file in files:
        parts = os.path.basename(file).split('_')
        tissue = '_'.join(parts[:2])
        df1 = pd.read_csv(os.path.join(directory, file), sep='\t', 
                          names=["Chr","Pos","Ref","Read Depth","Reads", "Quality"])
        
        A, C, G, T, Tis = [], [], [], [], []
        
        for reads, quality in zip(df1['Reads'], df1['Quality']): 
            a, c, g, t = 0, 0, 0, 0
            check_quality = len(reads) == len(quality)
            for idx, nucleotide in enumerate(reads):
                if nucleotide.upper() == 'A' and (not check_quality or quality[idx] not in low_quality_symbols):
                    a += 1
                elif nucleotide.upper() == 'C' and (not check_quality or quality[idx] not in low_quality_symbols):
                    c += 1
                elif nucleotide.upper() == 'G' and (not check_quality or quality[idx] not in low_quality_symbols):
                    g += 1
                elif nucleotide.upper() == 'T' and (not check_quality or quality[idx] not in low_quality_symbols):
                    t += 1
            
            A.append(a)
            C.append(c)
            G.append(g)
            T.append(t)
            Tis.append(tissue)

        df1['A'] = A
        df1['C'] = C
        df1['G'] = G
        df1['T'] = T
        df1['Tissue'] = Tis
        df_list.append(df1)  # Add the DataFrame to the list

    final_df = pd.concat(df_list, ignore_index=True)  # Concatenate all the DataFrames in the list
    
    return final_df

df = CountUp()
df.to_csv('./ASDEG/results/samtools/Nuc_countup')
