
import pandas as pd
import time

start2 = time.time()
similarity = 1
interv = -0.01
dfparame = pd.DataFrame()
activelist= []
retrievesimilar= []
Redundantretrieve= []
similaritylist= []
alllist= []
with open(r"in-silico toxicity and drug interaction prediction\integration\similar_comp_properties_sum.csv") as fin:
    df1 = pd.read_csv(fin)

while similarity <=1 and similarity >=0.01:
    print(df1.shape[0])
    RepreComList= list(set(df1["Active_compound_name"]))[0:]
    mask = (df1["Similarity_score"] >= similarity) & (df1["Active_compound_name"].isin(RepreComList))
    df2 = df1[mask]
    
    print(df2.shape[0])
    df2 = df2.drop(['Similarity_score'],axis= 1)
    df3 = df2.drop_duplicates(subset=None, keep='first', inplace=False)

    print("number of Actives "+ str(list(set(df2["Active_compound_name"]))))
    activelist.append(len(list(set(df2["Active_compound_name"]))))

    print("number of Retireves "+ str(df2.shape[0]))
    retrievesimilar.append(df2.shape[0])

    Redundantretrieve.append(round(((df2.shape[0]-df3.shape[0])/df2.shape[0]),3))
    print("Redundancy ration "+str(round(((df2.shape[0]-df3.shape[0])/df2.shape[0]),3)))
    similaritylist.append(round(similarity, 3))
    similarity += interv

dfparame= pd.DataFrame(list(zip(activelist, retrievesimilar, Redundantretrieve, similaritylist)), columns = ["number_of_Actives","number_of_Similar_compounds","Redundancy_rate","Similarity_threshold"])
dfparame.to_csv(r"in-silico toxicity and drug interaction prediction\integration\parameter research\parameter_similar_compound_properies.csv")
