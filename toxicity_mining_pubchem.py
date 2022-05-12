# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 15:13:48 2021

@author: User
"""
from requests_html import HTMLSession
import pandas as pd
import requests, re, os, argparse, time, json

def get_toxicity_information(cid):
    session = HTMLSession()
    new = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/'+str(cid)+'/JSON/?heading=Toxicological Information'
    response = session.get(new)
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
        print(time.time())
        df = pd.DataFrame([retrieve], columns= headlist)
        return df
    
if __name__ == "__main__":
            
    start = time.time()
    df1 = pd.DataFrame()
    with open() as fin:
        df = pd.read_csv(fin)
        com_cids = df['PubChem_CID_of_similar_compounds']
        com_name = df["Active_compound_name"]
        similarityscore = df["Similarity_score"]
    for cid, name, score in zip(com_cids[0:], com_name[0:], similarityscore[0:]):
        df2 = get_toxicity_information(cid)
        df2["PubChem_CID_of_similar_compounds"] = cid
        df2["Active_compound_name"] = name
        df2["Similarity_score"] = score
        df1 = df1.append(df2, ignore_index = True)
        
    df1.to_csv()
    end = time.time()
    print(df2.tail(5))
    print('execution time was {} minutes'.format(round((end-start)/60, 2)))