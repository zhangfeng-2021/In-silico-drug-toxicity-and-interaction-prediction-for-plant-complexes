from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException
from selenium.common.exceptions import UnexpectedAlertPresentException
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import Select
import pyautogui as pag

import time
import re
import pandas as pd

### DEFINE FUNCTIONs ### (Updated on 2022 07 08)

def SwissSimilarCrawler (smiles, CpdName): 
    option = Options()
    option.add_argument("--headless")
    SwissUrl = 'http://www.swisssimilarity.ch/' # URL of mining page (Updated on 2022 07 08)
    platform = 'SwissSimilarity'
    driver = webdriver.Chrome()
    try:
        driver.get(SwissUrl) # browsing to mining page
        time.sleep(2)
        SearchField = driver.find_element_by_name('smiles')  # locate smiles input fied
        time.sleep(0.5)
        SearchField.send_keys(smiles) # send smiles string to smiles input field
        time.sleep(0.5)
        driver.find_element_by_name('compoundClasses').click()
        time.sleep(0.5)
        # menu = Select(driver.find_elements_by_xpath("//select[@id='compoundClasses']"))
        # print("find menue")
        # menu.select_by_visible_text("Bioactive")
        # print("select bioactive")
        # Bioactive = driver.find_elements_by_xpath("//select/option[@value='2']")
        # print("find bioactie")
        # time.sleep(2)
        # print(pag.position())
        pag.moveTo(573,824) # The coordinate depends on your own situation
        pag.click(573,824)
        driver.find_element_by_xpath("//input[@value='CHEMBL_act___9']").click() # choose the libray and method
        time.sleep(0.5)
        SearchField.submit() # submit filled form
        time.sleep(0.5)
        # print(driver.window_handles)
        driver.switch_to_window(driver.window_handles[1])   
        df = pd.DataFrame() #vector of total similar compounds information of one active compound
        try:
            WebDriverWait(driver, 60).until(EC.visibility_of_all_elements_located((By.XPATH, "//tr"))) # wait for the result page to be loaded
            # print("table resrult visiable")
            CurrUrl = driver.current_url # get url of result page
            # print(CurrUrl)
            totals_rows = driver.find_elements_by_xpath("html/body/div[1]/div[2]/div[2]/table[1]/tbody[1]/tr")
            # print("totals_rows")
            listid = []
            listscore = []
            listsmiles = []
            _count = 0
            _a = 0
            for row in totals_rows[_a:]:
                rawsmilespool = row.get_attribute('innerHTML')
                # print(rawsmilespool)
                patterna = r"\<script\>(.+)\<\/script\>"
                rawsmiles = re.findall(patterna, rawsmilespool)
                # print(rawsmiles[0:1])
                _i = 0
                while _i < len(rawsmiles):
                    # print(rawsmiles[_i])
                    patternb = r"\=.*"
                    # print(re.findall(patternb, rawsmiles[_i]))
                    smiles = re.findall(patternb, rawsmiles[_i])[0][2:-2]
                    # print(smiles)
                    listsmiles.append(smiles)
                    _i += 1
                # print(row.text)
                row = row.text.split("\n")
                # print("row: "+ str(_count+_a+1), row)
                __i,__j = 0,1            
                while __i <= len(row)-1:
                    conta = row[__i]
                    listid.append(conta)
                    # print(conta)
                    __i += 2
                while __j <= len(row)-1:
                    contb = row[__j]
                    # print(contb)
                    listscore.append(contb[8:])
                    __j += 2
                _count += 1
                
            print("====================================")
            print(len(listid), len(listscore))
            print(listid[0:2], listscore[0:2], listsmiles[0:2])
            df = pd.DataFrame(zip(listid, listscore, listsmiles), columns = ["CHEMBL_Similar_Compound_ID", "Similarity_score", "SMILES_of_similar"])
            print("Data of one active compound stored")
            df["Active_compound_name"] = CpdName
            print("Active compound name insertion")
            df["Platform"] = platform
            print("Platform insertion")
            # dfKeep = df[df.columns[2]] > 0.09999 # remove all below probability of 0.1
            # df = df[dfKeep]
            return df
        except TimeoutException:
            CurrUrl = driver.current_url
            df = pd.DataFrame(columns=['CHEMBL_Similar_Compound_ID','Similarity_score','SMILES_of_similar'])
            df = df.append({'Active_compound_name':CpdName, # create a row with name; platform; url of result table
                'Platform': platform,
                "Issue_url": CurrUrl},
                ignore_index=True)
            return df
        except UnexpectedAlertPresentException:
            alert = driver.switch_to.alert
            df = pd.DataFrame(columns=['v','Similarity_score','SMILES_of_similar'])
            df = df.append({'Active_compound_name': CpdName, # creates row with name; platform; url of result table
                'Platform':platform, 
                "Issue": "UnexpectedAlertPresentException"},
                ignore_index=True)
            alert.accept()
            return df    
    except:
        raise
        
    finally:
        driver.close()
        time.sleep(0.5)

if __name__ == "__main__":        
    start = time.time()
    df1 = pd.read_csv(r"in-silico toxicity and drug interaction prediction\Property_active\Active_compounds_properties2022_07_08.csv")
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
    df2.to_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\similar_comp_pool_swiss_2022_07_08.csv") 
    seta = set(activenamelist)-set(retrinamelist) 
    df3.to_csv(r"in-silico toxicity and drug interaction prediction\Similar compounds mining\no_similar_retrieve_2022_07_08.csv")    
    end = time.time()
    print('Similar not found:', seta)
    print('execution time was {} minutes'.format(round((end-start)/60, 2)))
