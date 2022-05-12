
import pandas as pd
import time

start2 = time.time()
similarity = 1
interv = -0.01
interval = -0.11
dfparame = pd.DataFrame()
activelist= []
retrievesimilar = []
Redundantretrieve = []
similaritylist = []

# dftoxi = pd.DataFrame()
# while similarity <=1 and similarity >0.02:
    # with open(r"Toxi_all_infor_similar_"+str((round((similarity+interv),3), round(similarity,3)))+".csv", 
    #           encoding = None, errors="ignore") as fin:
#         dfin = pd.read_csv(fin)
#     print(dfin.shape)
#     dftoxi = dftoxi.append(dfin)
#     similarity += interval
# dftoxi= dftoxi.sort_values(by= ['Active_compound_name', 'Similarity_score'], ascending= [False, False])
# dftoxi.to_csv(r"Toxi_infor_sum.csv")
with open(r"Toxi_infor_sum.csv", encoding = None, errors="ignore") as fin:
    df1 = pd.read_csv(fin)

while similarity <= 1 and similarity >= 0.01:
    print(df1.shape[0])
    RepreComList = list(set(df1["Active_compound_name"]))[0:]
    mask = (df1["Similarity_score"] >= similarity) & (df1["Active_compound_name"].isin(RepreComList))
    df2 = df1[mask]

    print(df2.shape[0])
    df2 = df2.drop(['Similarity_score'],axis = 1)
    df3 = df2.drop_duplicates(subset=None, keep ='first', inplace=False)

    print("number of Actives "+ str(list(set(df2["Active_compound_name"]))))
    activelist.append(len(list(set(df2["Active_compound_name"]))))

    print("number of Retireves "+ str(df2.shape[0]))
    retrievesimilar.append(df2.shape[0])
    if df2.shape[0] == 0:
        Redundantretrieve.append(None)
    else:
        Redundantretrieve.append(round(((df2.shape[0]-df3.shape[0])/df2.shape[0]),3))
        print("Redundancy ratio "+str(round(((df2.shape[0]-df3.shape[0])/df2.shape[0]),3)))
    similaritylist.append(round(similarity, 3))
    similarity += interv

dfparame= pd.DataFrame(list(zip(activelist, retrievesimilar, Redundantretrieve, similaritylist)), 
                       columns = ["number_of_Actives","Retrieves_of_Toxi_information",
                                  "Redundancy_rate_of_toxi_mining","Similarity_threshold"])
dfparame.to_csv(r"parameter_toxi_infor.csv")
