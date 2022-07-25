# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 14:59:10 2021

@author: User
"""
import pandas as pd
import time, math
from requests_html import HTMLSession
from lxml.html import fromstring
import requests, re, os, argparse,json
import urllib.parse
import urllib.request
import numpy as np
import pyautogui as pag
import random

## toxicity mining function definition ##   
def get_toxicity_information(cid):
    session = HTMLSession()
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/'+str(int(cid))+'/JSON/?heading=Toxicological Information'
    response = session.get(url)
    tox = response.html.html
    tox = json.loads(tox)
    #print(tox)
    headlist = []
    retrieve = []
    df = pd.DataFrame()
    
    try:
        i = 0
        while tox["Record"]["Section"][0]["Section"][0]["Section"][i]["TOCHeading"]:
            head= tox["Record"]["Section"][0]["Section"][0]["Section"][i]["TOCHeading"]
            headlist.append(head)
            try:
                j = 0
                content = []
                while tox["Record"]["Section"][0]["Section"][0]["Section"][i]["Information"][j]:
                    strings = tox["Record"]["Section"][0]["Section"][0]["Section"][i]["Information"][j]["Value"]["StringWithMarkup"][0]["String"]
                    content.append(strings) 
                    j += 1
            except:
                retrieve.append(content)
            i += 1
    except:
        print(len(retrieve), headlist)
        df = pd.DataFrame([retrieve], columns= headlist)
        return df

##  ##    
if __name__ == "__main__":
    start = time.time() 
    similarity = 1
    interv = -0.11
    while similarity <=1 and similarity >=0.1:
        with open(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\Similar compounds property\similar_all_properties_"+str((round((similarity+interv),3), round(similarity,3)))+".csv", encoding = None, errors="ignore") as fin:
            dfsimilar = pd.read_csv(fin) # Read similar compound dataframe
        RepreComList = list(set(dfsimilar["Active_compound_name"]))[0:]+[""] # Set interested compounds
        #print(RepreComList)
        mask = ((dfsimilar["Similarity_score"] > (similarity + interv)) & (dfsimilar["Similarity_score"] <= similarity) & (dfsimilar["Active_compound_name"].isin(RepreComList))) # Conditions of interested compound sorting
        dfsimilar = dfsimilar[mask] # Store the interested similar compounds in a new DataFrame file named dfsimilar
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
                print("Working on the compound "+ name + " out of "+ str(dfsimilar.shape[0]+1)+ " similar compounds")
                print("***************************")
                dfins = dfinhouse[dfinhouse["PubChem_CID"] == cid] # Check if it has been stored in the in-house vector
                print("PubChem_CID: "+ str(int(cid)))
                if cid in dfins["PubChem_CID"].values.tolist():# if it has been stored in the in-house vector
                    print("it has been stored in the in-house vector")
                    dfins = dfins.drop(labels = ["Active_compound_name","Similarity_score"], axis = 1)
                    print("Drop two coloumns which are not useful")
                    dfinst = dfins.drop_duplicates(subset= "PubChem_CID", keep='first')
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
                        print(i, cid)
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
