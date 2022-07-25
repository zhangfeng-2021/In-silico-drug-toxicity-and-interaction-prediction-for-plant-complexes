# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 15:25:12 2022

@author: User
"""
import pandas as pd
import pubchempy as pcp
import time
import pyautogui as pag
from similar_comp_crawler import SwissSimilarCrawler
import math
from toxicity_mining_pubchem import get_toxicity_information
import ast

start = time.time()

## Active compound property collection ##
with open(r"in-silico toxicity and drug interaction prediction\Active compounds_TCMSP\Active_compounds_pool.csv") as fin:
    df1 = pd.read_csv(fin)
compoundname= list(set(df1['Molecule Name']))[:]
df2 = pd.DataFrame()
try: 
    for com in compoundname:
        compounds = pcp.get_compounds(com,'name')
        for coms in compounds:
                df3 = pcp.compounds_to_frame(coms)
                df3['Molecule Name'] = com
                df2 = df2.append(df3)
                print(com,coms.isomeric_smiles)
except:
    print(com,'no hit')
    df2 = df2.append([com, 'no hit'])
df2.rename(columns = {'cid':'PubChem_CID', 'Molecule Name':'Active_compound_name'}, inplace = True)
retrinamelist = df2['Active_compound_name']
df2.to_csv(r"in-silico toxicity and drug interaction prediction\Property_active\Active_compounds_properties.csv") 
seta = set(compoundname)-set(retrinamelist) 
df4 = pd.DataFrame(list(seta), columns=['Compound not found in Pubchem'])    
df4.to_csv(r"in-silico toxicity and drug interaction prediction\Property_active\Active_compounds_properties_no_retrie.csv")    
end = time.time()
print('Compound not found in Pubchem:', seta)
print('execution time was {} minutes'.format(round((end-start)/60, 2)))

## Similar compounds mining ##
df1 = pd.read_csv(r"in-silico toxicity and drug interaction prediction\Property_active\Active_compounds_properties.csv")
compoundsmiles = df1['isomeric_smiles'][:]
activenamelist = df1["Active_compound_name"][:]
df2 = pd.DataFrame()
df3 = pd.DataFrame(columns =['no_similar_compound_found',"isomeric_smiles"])
i = 0
for smiles, CpdName in zip(compoundsmiles[:], activenamelist[:]):
    try:
        retrieve = SwissSimilarCrawler (smiles, CpdName)
        df2 = df2.append(retrieve)
        print(time.time())
    except:
        print(CpdName,smiles,'no hit', time.time())
        df3.loc[i] = ([CpdName, smiles])
        i += 1
print(df2.head(2))                   
retrinamelist = df2['Active_compound_name']
df2.to_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\similar_comp_pool.csv") 
seta = set(activenamelist)-set(retrinamelist) 
df3.to_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\no_similar_retrieve.csv")    
end = time.time()
print('Similar not found:', seta)
print('execution time was {} minutes'.format(round((end-start)/60, 2)))


### Similar compound property and drug toxicity interaction mining ####

similarity = 1
interv = -0.11
start1 = time.time()
while similarity <=1 and similarity >=0.10:
    
# Similar compound property mining #
    df1 = pd.read_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\similar_comp_pool.csv")
    RepreComList = list(set(df1["Active_compound_name"]))[0:]+[""]
    # print(RepreComList)
    mask = ((df1["Similarity_score"] > (similarity+interv)) & (df1["Similarity_score"] <= similarity) & (df1["Active_compound_name"].isin(RepreComList)))
    df1 = df1[mask]
    # print(df1)
    compoundsmiles = df1['SMILES_of_similar']
    activenamelist = df1["Active_compound_name"]
    CHEMBLlist = df1["CHEMBL_Similar_Compound_ID"]
    Similarity_scorelist = df1["Similarity_score"]
    dfstore = pd.DataFrame(columns =['Active_compound_name', "CHEMBL_Similar_Compound_ID", "PubChem_CID", "Similarity_score"])
    dfno = pd.DataFrame(columns =['Similar not found in Pubchem',"SMILES_of_similar","CHEMBLID"])
    dfinhouse = pd.read_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\in house vector\similar_comp_total_in_house.csv", index_col= 0)
    i = 0
    try:
        for smiles, name, CHEMBLID, prob in zip(compoundsmiles[0:], activenamelist[0:], CHEMBLlist[0:], Similarity_scorelist[0:]):
            a,b = pag.position()
            time2 = time.time()
            timecost1 = time2 - start1
            print(math.ceil(timecost1/30))
            if math.ceil(timecost1/30) %100 == 0:
                time.sleep(20)
            else:
                print("\nSimilar compound property mining")
                print("=====================================")
                # print(df1.shape)
                print("Working on the compound "+ str(i+1)+": "+str(CHEMBLID)+ ", out of "+ str(df1.shape[0])+ " similar compounds")
                # print(CHEMBLID)
                print("Similarity range: "+str((round((similarity+interv),3), round(similarity,3))))
            
            dfins = dfinhouse[dfinhouse["CHEMBL_Similar_Compound_ID"]== CHEMBLID]
            if CHEMBLID in dfinhouse["CHEMBL_Similar_Compound_ID"].values.tolist():
                print("This compound has been stored in former vector")
                # print(dfins)
                dfins = dfins.drop(labels = ["Active_compound_name","Similarity_score"], axis= 1)
                print("Drop two coloumns which are not useful")
                # print(dfins)
                dfins = dfins.drop_duplicates(subset= "CHEMBL_Similar_Compound_ID", keep='last')
                print("Drop duplicates")
                # print(dfins)
                dfins["Active_compound_name"] = name
                dfins["Similarity_score"] = prob
                dfstore = dfstore.append(dfins, ignore_index = True)
                print("It has been stored successfully")
                i += 1
            else:
                try:
                    cidlist = pcp.get_cids(CHEMBLID,'name')
                    if cidlist and cidlist[0] != 0 :
                        # print("could by name")
                        # print(cidlist)
                        for cids in cidlist:
                            print("PubChem_CID: " + str(cids))
                            df = pd.DataFrame([[name, CHEMBLID, cids, prob]],columns =['Active_compound_name', "CHEMBL_Similar_Compound_ID", "PubChem_CID", "Similarity_score"])
                            dfstore = dfstore.append(df, ignore_index = True)
                            print(time.time(), "Retrieve success by name")
                            i += 1
                    else:
                        cidlist = pcp.get_cids(smiles,'smiles')
                        if cidlist:
                            for cids in cidlist:
                                df = pd.DataFrame([[name, CHEMBLID, cids, prob]],columns =['Active_compound_name', "CHEMBL_Similar_Compound_ID", "PubChem_CID", "Similarity_score"])            
                                # print(str(cids))
                                dfstore = dfstore.append(df, ignore_index = True)
                                print(time.time(),"Retrieve success by ID-",CHEMBLID)
                                i += 1
                        else:
                            print(name,smiles,'no hit', time.time())
                            dfno.loc[i] = ([name, smiles, CHEMBLID])
                            i += 1
                except:
                    print(name,smiles,'no hit', time.time())
                    dfno.loc[i] = ([name, smiles, CHEMBLID])
                    i += 1
        dfinhouse = pd.read_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\in house vector\similar_comp_total_in_house.csv", index_col= 0) # Access the in-house vector
        dfinhouse = dfinhouse.append(dfstore, ignore_index = True) 
        dfinhouse = dfinhouse.drop_duplicates(subset = ["Active_compound_name","PubChem_CID"])
        dfinhouse.to_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\in house vector\similar_comp_total_in_house.csv")
        retrinamelist = dfstore['Active_compound_name']
        dfstore.to_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\Similar compounds property\similar_all_properties_"+str((round((similarity+interv),3), round(similarity,3)))+".csv") 
        seta = set(activenamelist)-set(retrinamelist) 
        dfno.to_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\Similar compounds property\no retrieve\similar_compound_pubchem_not_found_"+str(round(similarity, 3))+".csv")    
        end1 = time.time()
        # print('Compound not found in Pubchem:', seta)
        print('execution time was {} minutes'.format(round((end1-start1)/60, 2)))
    except:        
        dfinhouse = pd.read_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\in house vector\similar_comp_total_in_house.csv", index_col= 0) # Access the in-house vector
        dfinhouse = dfinhouse.append(dfstore, ignore_index = True) 
        dfinhouse = dfinhouse.drop_duplicates(subset = ["Active_compound_name","PubChem_CID"])
        dfinhouse.to_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\in house vector\similar_comp_total_in_house.csv")
        retrinamelist = dfstore['Active_compound_name']
        dfstore.to_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\Similar compounds property\similar_all_properties_"+str((round((similarity+interv),3), round(similarity,3)))+".csv") 
        seta = set(activenamelist)-set(retrinamelist) 
        dfno.to_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\Similar compounds property\no retrieve\similar_compound_pubchem_not_found_"+str(round(similarity, 3))+".csv")    
        end1 = time.time()
        # print('Compound not found in Pubchem:', seta)
        print("Completed")
        print('execution time was {} minutes'.format(round((end1-start1)/60, 2)))
        raise
        
# Drug toxicity interaction mining #
    with open(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\Similar compounds property\similar_all_properties_"+str((round((similarity+interv),3), round(similarity,3)))+".csv", encoding = None, errors="ignore") as fin:
        dfsimilar = pd.read_csv(fin) # Read similar compound dataframe
    RepreComList = list(set(dfsimilar["Active_compound_name"]))[0:]+[""] # Set interested compounds
    #print(RepreComList)
    mask = (dfsimilar["Similarity_score"] > (similarity + interv)) & (dfsimilar["Similarity_score"] <= similarity) & (dfsimilar["Active_compound_name"].isin(RepreComList)) # Conditions of interested compound sorting
    dfsimilar = dfsimilar[(mask)] # Store the interested similar compounds in a new DataFrame file named dfsimilar
    compoundsmiles = dfsimilar['isomeric_smiles']
    CHEMBLlist = dfsimilar["CHEMBL_Similar_Compound_ID"]
    com_cids = dfsimilar["PubChem_CID"]
    com_name = dfsimilar["Active_compound_name"]
    similarityscore = dfsimilar["Similarity_score"]
    dfstore = pd.DataFrame()
    dfinhouse = pd.read_csv(r"in-silico toxicity and drug interaction prediction\Toxicity and interaction mining\in house vector\Toxi_infor_in_house.csv", index_col= 0) # Access the in-house vector
    print("Find in-house toxicity information vector", dfinhouse.shape)
    
    try:
        i = 0
        for cid, name, similar in zip(com_cids[0:], com_name[0:], similarityscore[0:]): # Iterate every interested similar compound
            print("\nDrug toxicity interaction mining")
            print("Working on the compound "+ name+ " out of "+ str(dfsimilar.shape[0]+1)+ " similar compounds")
            print("***************************")
            dfins = dfinhouse[dfinhouse["PubChem_CID"] == cid] # Check if it has been stored in the in-house vector
            print("PubChem_CID: "+ str(int(cid)))
            if cid in dfins["PubChem_CID"].values.tolist():# if it has been stored in the in-house vector
                print("it has been stored in the in-house vector")
                dfins = dfins.drop(labels = ["Active_compound_name","Similarity_score"], axis = 1)
                print("Drop two coloumns which are not useful")
                dfinst = dfins.drop_duplicates(subset= "PubChem_CID", keep='last')
                print("Drop duplicates")
                dfinst["Active_compound_name"] = name
                dfinst["Similarity_score"] = similar
                dfstore = dfstore.append(dfinst, ignore_index = True)
                print("Similarity range: "+str((round((similarity+interv),3), round(similarity,3))))
                print("It has been appended")
                i += 1
            else: # if it has not been stored in the in-house vector
                c,d = pag.position()
                df1 = get_toxicity_information(cid)
                print("Search toxicity information")
                if df1.empty:
                    print(str(cid)+" "+str(name)+" no retrieve")
                    print("Similarity range: "+str((round((similarity+interv),3), round(similarity,3))))
                    # df1["PubChem_CID"] = cid
                    # df1["Active_compound_name"] = name
                    # df1["Similarity_score"] = similar
                    # dfstore = dfstore.append(df1, ignore_index = True)
                    # print(dfstore.shape)
                    print("No toxicity and interaction information of "+ str(i) +" similar compounds have been collected")
                else:
                    df1["PubChem_CID"] = cid
                    print(cid)
                    df1["Active_compound_name"] = name
                    df1["Similarity_score"] = similar
                    dfstore = dfstore.append(df1, ignore_index = True)
                    print(dfstore.shape)
                    print("Similarity range: "+str((round((similarity+interv),3), round(similarity,3))))
                    print("Toxicity and interaction information of "+ str(i) +" similar compounds have been collected")
                presenttime = time.time()
                timecost1 = presenttime - start
                print("Time cost: "+str(round(timecost1, 1))+" s")
                print(math.ceil(timecost1/10))
                if math.ceil(timecost1/10) %50 == 0:
                    time.sleep(15)
                else:
                    print("Missions of "+ str(i) +" similar compounds have been completed")
                i += 1
                print("******************************************") 
        dfstore.to_csv(r"in-silico toxicity and drug interaction prediction\Toxicity and interaction mining\Toxi_all_infor_similar_"+str((round((similarity+interv),3), round(similarity,3)))+".csv", encoding = None, errors="ignore")                    
        dfinhouse = pd.read_csv(r"in-silico toxicity and drug interaction prediction\Toxicity and interaction mining\in house vector\Toxi_infor_in_house.csv", index_col= 0) # Access the in-house vector
        dfinhouse = dfinhouse.append(dfstore, ignore_index = True) 
        dfinhouse = dfinhouse.drop_duplicates(subset = ["Active_compound_name","PubChem_CID","Similarity_score"], keep='last')
        dfinhouse.to_csv(r"in-silico toxicity and drug interaction prediction\Toxicity and interaction mining\in house vector\Toxi_infor_in_house.csv")
        end2 = time.time()
        print('execution time was {} minutes'.format(round((end2-start)/60, 2)))
    except:
        dfstore.to_csv(r"in-silico toxicity and drug interaction prediction\Toxicity and interaction mining\Toxi_all_infor_similar_"+str((round((similarity+interv),3), round(similarity,3)))+".csv", encoding = None, errors="ignore")            
        dfinhouse = pd.read_csv(r"in-silico toxicity and drug interaction prediction\Toxicity and interaction mining\in house vector\Toxi_infor_in_house.csv", index_col= 0) # Access the in-house vector
        dfinhouse = dfinhouse.append(dfstore, ignore_index = True)
        dfinhouse = dfinhouse.drop_duplicates(subset = ["Active_compound_name","PubChem_CID","Similarity_score"], keep='last')
        dfinhouse.to_csv(r"in-silico toxicity and drug interaction prediction\Toxicity and interaction mining\in house vector\Toxi_infor_in_house.csv")
        end2 = time.time()
        print('execution time was {} minutes'.format(round((end2-start)/60, 2)))
        raise
    similarity += interv
print("Toxicity raw data collected")

## Data integration ##
start = time.time()
dftoxitotal = pd.DataFrame()
dfsimitotal = pd.DataFrame()
start = time.time() 
similarity = 1
interv = -0.11
while similarity <=1 and similarity >=0.10:
    print("Similarity range: "+str((round((similarity+interv),3), round(similarity,3))))
    with open(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\Similar compounds property\similar_all_properties_"+str((round((similarity+interv),3), round(similarity,3)))+".csv", encoding = None, errors="ignore") as fin:
        dfsimilar = pd.read_csv(fin, index_col= 0) # Read similar compound dataframe
    with open(r"in-silico toxicity and drug interaction prediction\Toxicity and interaction mining\Toxi_all_infor_similar_"+str((round((similarity+interv),3), round(similarity,3)))+".csv", encoding = None, errors="ignore") as fin:
        dftoxi = pd.read_csv(fin, index_col= 0) # Read similar compound dataframe
    dfinhouse = pd.read_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\in house vector\similar_comp_total_in_house.csv", index_col= 0)
    dfinhouse = dfinhouse.append(dfsimilar, ignore_index = True)
    dfinhouse = dfinhouse.drop_duplicates(subset = ["Active_compound_name","PubChem_CID"])
    dfinhouse.to_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\in house vector\similar_comp_total_in_house.csv") 
    print("More similar compounds properties added to in-house vector")
    dfinhouse = pd.read_csv(r"in-silico toxicity and drug interaction prediction\Toxicity and interaction mining\in house vector\Toxi_infor_in_house.csv", index_col= 0)
    dfinhouse = dfinhouse.append(dftoxi, ignore_index = True)
    dfinhouse = dfinhouse.drop_duplicates(subset = ["Active_compound_name","PubChem_CID","Similarity_score"], keep='last')    
    dfinhouse.to_csv(r"in-silico toxicity and drug interaction prediction\Toxicity and interaction mining\in house vector\Toxi_infor_in_house.csv")
    print("More toxicity and drug interaction information added to in-house vector")
    dftoxitotal = dftoxitotal.append(dftoxi, ignore_index = True)
    dfsimitotal = dfsimitotal.append(dfsimilar, ignore_index = True)
    similarity += interv 
dftoxitotal = dftoxitotal.sort_values(['Active_compound_name', 'Similarity_score'], ascending=[True, False], axis = 0)
dfsimitotal = dfsimitotal.sort_values(['Active_compound_name', 'Similarity_score'], ascending=[True, False], axis = 0)
dftoxitotal.to_csv(r"in-silico toxicity and drug interaction prediction\integration\Toxi_infor_sum.csv")
mask = (dftoxitotal["Drug Induced Liver Injury"].notna()|dftoxitotal["Hepatotoxicity"].notna() |dftoxitotal["Evidence for Carcinogenicity"].notna() |dftoxitotal["Carcinogen Classification"].notna())
dftoxitotal = dftoxitotal[mask]
dftoxitotal.to_csv(r"in-silico toxicity and drug interaction prediction\integration\Toxi_infor_sum_similarity_threshold_analysis_data.csv")
dfsimitotal.to_csv(r"in-silico toxicity and drug interaction prediction\integration\similar_comp_properties_sum.csv")
end = time.time()
print('execution time was {} minutes'.format(round((end-start)/60, 2)))

## Threshold analysis ## 

# Drug interaction information split #
threshold = 0.6171
dfintera = pd.DataFrame(columns = ["Active_compound_name","PubChem CID of similar compounds","Similarity_score",
                              "Interactions","Drug_interacted", "Activities","Roles"])
dfins = pd.DataFrame(columns = ["Active_compound_name","PubChem CID of similar compounds","Similarity_score",
                              "Interactions","Drug_interacted", "Activities","Roles"])

with open(r"in-silico toxicity and drug interaction prediction\integration\Toxi_infor_sum.csv", encoding = None, errors="ignore") as fin:
    df1 = pd.read_csv(fin)
container = []
df1 = df1[df1["Interactions"].notna()]
mask = df1["Similarity_score"] > threshold
df1 = df1[mask] 
for inter, cid, name, score in zip(df1["Interactions"], df1["PubChem_CID"], 
                                   df1["Active_compound_name"], df1["Similarity_score"]):
    inter = ast.literal_eval(inter)                                                            
    for j in inter:
        print(j)
        container.append(name)
        container.append(cid)
        container.append(score)
        container.append(j)
        dfins = pd.DataFrame([container], columns = ["Active_compound_name","PubChem CID of similar compounds",
                                                   "Similarity_score","Interactions"])
        dfintera = dfintera.append(dfins, ignore_index= True)
        container = []
dfintera.drop_duplicates(subset=(["Active_compound_name", "PubChem CID of similar compounds","Interactions"]), keep= 'last')
dfintera.to_csv(r"in-silico toxicity and drug interaction prediction\Toxicity and interaction mining\Drug_interaction\drug_interaction_"+str(threshold)+".csv")
print("Interaction data splitted")

## Parameter analysis##

#Similar compound related#
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

#Toxicity related#
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
with open(r"in-silico toxicity and drug interaction prediction\integration\Toxi_infor_sum.csv", encoding = None, errors="ignore") as fin:
    df1 = pd.read_csv(fin)

while similarity <= 1 and similarity >= 0.01:
    print(df1.shape[0])
    RepreComList = list(set(df1["Active_compound_name"]))[0:]
    mask = (df1["Similarity_score"] >= similarity) & (df1["Active_compound_name"].isin(RepreComList))
    df2 = df1[mask]

    print(df2.shape[0])
    df2 = df2.drop(['Similarity_score'], axis = 1)
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
dfparame.to_csv(r"in-silico toxicity and drug interaction prediction\integration\parameter research\parameter_toxi_infor.csv")

#Active compounds predicted based on similarity threshold#
with open(r"in-silico toxicity and drug interaction prediction\integration\Toxi_infor_sum.csv", encoding = None, errors="ignore") as fin:
    df1 = pd.read_csv(fin)
    
similarity = 1
RepreComList1= list(set(df1["Active_compound_name"]))[0:]
mask1 = (df1["Similarity_score"] >= similarity) & (df1["Active_compound_name"].isin(RepreComList1))
df2 = df1[mask1]
activenumber1 = len(set(df2["Active_compound_name"]))
print("Number of active compounds with toxicity retrieved: ", str(activenumber1), " with the similarity threshold of "+ str(similarity))
    

similarity = 0.759 
RepreComList1= list(set(df1["Active_compound_name"]))[0:]
mask1 = (df1["Similarity_score"] > similarity) & (df1["Active_compound_name"].isin(RepreComList1))
df2 = df1[mask1]
activenumber1 = len(set(df2["Active_compound_name"]))
print("Number of active compounds with toxicity retrieved: ", str(activenumber1), " with the similarity threshold of "+ str(similarity))

similarity = 0.6171 
RepreComList1= list(set(df1["Active_compound_name"]))[0:]
mask1 = (df1["Similarity_score"] > similarity) & (df1["Active_compound_name"].isin(RepreComList1))
df2 = df1[mask1]
activenumber1 = len(set(df2["Active_compound_name"]))
df3 = df2.drop_duplicates(subset=["Active_compound_name", "PubChem_CID"], keep='first', inplace=False)
df3.to_csv(r"in-silico toxicity and drug interaction prediction\integration\prediction_"+str(similarity)+".csv")
print("Number of active compounds with toxicity retrieved: ", str(activenumber1), " with the similarity threshold of "+ str(similarity))


with open(r"in-silico toxicity and drug interaction prediction\integration\similar_comp_properties_sum.csv", encoding = None, errors="ignore") as fin:
    df1 = pd.read_csv(fin)
    
similarity = 1 
RepreComList1= list(set(df1["Active_compound_name"]))[0:]
mask1 = (df1["Similarity_score"] >= similarity) & (df1["Active_compound_name"].isin(RepreComList1))
df2 = df1[mask1]
activenumber1 = len(set(df2["Active_compound_name"]))
print("Number of active compounds with similar compounds retrieved: ", str(activenumber1), " with the similarity threshold of "+ str(similarity))

similarity = 0.759 
RepreComList1= list(set(df1["Active_compound_name"]))[0:]
mask1 = (df1["Similarity_score"] > similarity) & (df1["Active_compound_name"].isin(RepreComList1))
df2 = df1[mask1]
activenumber1 = len(set(df2["Active_compound_name"]))
print("Number of active compounds with similar compounds retrieved: ", str(activenumber1), " with the similarity threshold of "+ str(similarity))

similarity = 0.6171 
RepreComList1= list(set(df1["Active_compound_name"]))[0:]
mask1 = (df1["Similarity_score"] > similarity) & (df1["Active_compound_name"].isin(RepreComList1))
df2 = df1[mask1]
activenumber1 = len(set(df2["Active_compound_name"]))
print("Number of active compounds with similar compounds retrieved: ", str(activenumber1), " with the similarity threshold of "+ str(similarity))
