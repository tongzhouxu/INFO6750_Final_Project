from Bio import SeqIO
import pandas as pd
import random
import seaborn as sns

def calculate_kmer(size,sequence,seqname):
    kmer_dict = {}
    for i in range(len(sequence)-size+1):
        kmer_dict[sequence[i:i+size]] = 1
    return pd.DataFrame(kmer_dict,index=[seqname])

def similarity(center,target):
    sim_score = 0
    # higher score if two series have more same kmers
    for i in range(len(center)):
        if target[i] == center[i]:
            sim_score += 1/len(center)
    return sim_score

def kmodes(k,df,iteration):
    center_df = []
    # Initialization, randomly choose three unique rows from the dataset
    random_int_list = []
    while len(center_df) < k:
        random_int = random.randint(0, len(df.index)-1)
        if random_int not in random_int_list:
            center_df.append(df.iloc[random_int])

    group_index = [-1]*len(df.index) # a list of index# Iterations
    while iteration > 0:
        for i in range(len(df.index)):
            high_score = -1
            for j in range(k):
                score = similarity(center_df[j], df.iloc[i])
                if score > high_score:
                    high_score = score
                    group_index[i] = j
        
        # calculate the new center, i.e. the modes for each group
        for i in range(k):
            to_be_merged = []
            for j in range(len(group_index)):
                if group_index[j] == i:
                    to_be_merged.append(df.iloc[j])
            merged = pd.concat(to_be_merged,axis = 1)
            center_df[i] = merged.mode(axis = 1).iloc[:,0].astype(int)
        iteration = iteration -1
        print(group_index)
    return group_index

def predict(df,serotypes):
    predict_result = [-1]*len(df.index)
    count = 0
    for i in range(len(df.index)):
        high_score = -1
        for j in range(len(df.index)):
            if i == j:
                continue
            else:
                score = similarity(df.iloc[i], df.iloc[j])
                if score > high_score:
                    high_score = score
                    predict_result[i] = serotypes[j]
        if predict_result[i] == serotypes[i]:
            count += 1
    print("Predicted serotypes:", predict_result)
    print("Accuracy:", count/len(df.index)*100,'%')

def main():
    # Data Processing 1: identify all the kmers
    f = open('list.txt','r')
    name_list = f.readlines()
    file_list = []
    for i in name_list:
        file_list.append(i.strip('\n'))
    kmer_df_list = []
    for f in file_list:
        for record in SeqIO.parse(f,"fasta"):
            strain_name = record.id
            kmer_df = calculate_kmer(10,str(record.seq),record.id)
            kmer_df_list.append(kmer_df)
    df_merged = pd.concat(kmer_df_list).fillna(0).astype(int)
    df_merged.to_csv (r'kmer_df.csv', index = True, header=True)
    
    # Data Processing 2: filtering
    df = pd.read_csv('kmer_df.csv',header=0,index_col=0)
    nunique = df.apply(pd.Series.nunique)
    cols_to_drop = nunique[nunique == 1].index
    df = df.drop(cols_to_drop, axis=1)
    df.to_csv (r'kmer_df_filtered.csv', index = True, header=True)

    # Clustering
    df_filtered = pd.read_csv('kmer_df_filtered.csv',header=0,index_col=0)
    serotypes = ["Typhimurium" for i in range(10)] + ["Newport" for i in range(10)] + ["Enteritidis" for i in range(10)]
    group = kmodes(3,df_filtered,5)
    df_clustering_graph = pd.DataFrame([group,serotypes]).transpose()
    df_clustering_graph.columns = ['Group','Serotype']
    ax = sns.countplot(x='Group',hue='Serotype',data=df_clustering_graph)
    ax.figure.savefig('output')

    # Classification
    predict(df_filtered, serotypes)
    
main()