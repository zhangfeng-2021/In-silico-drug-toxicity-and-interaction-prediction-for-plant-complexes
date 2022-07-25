# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 11:40:41 2022

@author: User
"""
##  FC-SS & FE-SS annotation ##
import pandas as pd
import time
import re
import numpy as np

# string = "Overall evaluation: Cholesterol is not classifiable as to its carcinogenicity to humans (Group 3).']"
# string1 = "OVERALL EVALUATION: Group 3: The agent is not classifiable as to its carcinogenicity to humans.']"
# string2 = "['A4: Not classifiable as a human carcinogen.', 'TLV-A4']"
# string3 = " is reasonably anticipated to be a human carcinogen based on sufficient evidence of malignant tumor formation in multiple species of experimental animals.']"
# string4 = ""
# patternDILI = r", '(.*)-DILI-Concern"
# c = re.search(patternDILI, string)
# print(c[0])


##  Funtion definition ##
def DILI_annotation(df):
    dfassess = df[["Active_compound_name", 'Similarity_score',"PubChem_CID", "Drug Induced Liver Injury",
              "Hepatotoxicity", "Evidence for Carcinogenicity", "Carcinogen Classification"]]
    patternnoDILI = r"No-DILI-Concern"
    patternDILI = r"Most-DILI-Concern"
    patternDILInotsure1 = r"Less-DILI-Concern"
    patternDILInotsure2 = r"Ambiguous.DILI-concern"
    print(dfassess.head())
    for i in  dfassess["Drug Induced Liver Injury"]:
        if i == "":
            return np.NaN
        elif re.search(patternnoDILI, i, re.IGNORECASE):
            return -1
        elif re.search(patternDILInotsure1, i, re.IGNORECASE) or re.search(patternDILInotsure2, i, re.IGNORECASE):
            return 0
        elif re.search(patternDILI, i, re.IGNORECASE):
            return 1
        else:
            return "something wrong, please check mannually"
    
def Carcinogen_evidence_annotation (df):
    dfassess = df[["Active_compound_name", 'Similarity_score',"PubChem_CID", "Drug Induced Liver Injury",
              "Hepatotoxicity", "Evidence for Carcinogenicity", "Carcinogen Classification"]]
    print(dfassess.head())
    patternnotcarcino1 = r"Not classifiable"
    patternnotcarcino2 = r"Group 3"
    patterncarcinonotsure = r"Group 2B"
    patterniscarcino1 = r"Group 2A"
    patterniscarcino2 = r"reasonably anticipated to be"
    
    for i in dfassess["Evidence for Carcinogenicity"]:
        if i == "":
            return np.NaN
        elif re.search(patterncarcinonotsure, i, re.IGNORECASE):
            return 0
        elif re.search(patternnotcarcino2, i, re.IGNORECASE) or re.search(patternnotcarcino1, i, re.IGNORECASE):
            return -1
        elif re.search(patterniscarcino1, i, re.IGNORECASE) or re.search(patterniscarcino2, i, re.IGNORECASE):
            return 1
        else:
            return "something wrong, please check mannually"
    
def Carcinogen_classification_annotation (df):
    dfassess = df[["Active_compound_name", 'Similarity_score',"PubChem_CID", "Drug Induced Liver Injury",
              "Hepatotoxicity", "Evidence for Carcinogenicity", "Carcinogen Classification"]]
    print(dfassess.head())
    patternnotcarcino1 = r"Not classifiable"
    patternnotcarcino2 = r"Group 3"
    patterncarcinonotsure = r"Group 2B"
    patterniscarcino = r"Group 2A"    
    if dfassess["Carcinogen Classification"] == "":
        return np.NaN
    elif re.search(patterncarcinonotsure, dfassess["Carcinogen Classification"], re.IGNORECASE):
        return 0
    elif re.search(patternnotcarcino2, dfassess["Carcinogen Classification"], re.IGNORECASE) or re.search(patternnotcarcino1, 
                                                                                                          dfassess["Carcinogen Classification"], re.IGNORECASE):
        return -1
    elif re.search(patterniscarcino, dfassess["Carcinogen Classification"], re.IGNORECASE):
        return 1
    else:
        return "something wrong, please check mannually"
    
       
if __name__ == "__main__":
    dftoxianno = pd.read_csv(r"in-silico toxicity and drug interaction prediction\integration\Toxi_infor_sum_similarity_threshold_analysis_data.csv",index_col= 0)
    dftoxianno = dftoxianno[["Active_compound_name", 'Similarity_score',"PubChem_CID", "Drug Induced Liver Injury",
              "Hepatotoxicity", "Evidence for Carcinogenicity", "Carcinogen Classification"]]
    dftoxianno["Drug Induced Liver Injury"] = dftoxianno["Drug Induced Liver Injury"].replace()
    dftoxianno.loc[:,"DILI_annotation"] = dftoxianno.apply(DILI_annotation, axis = 1)
    dftoxianno.loc[:,"Carcinogen_evidence_annotation"] = dftoxianno.apply(Carcinogen_evidence_annotation, axis = 1)
    dftoxianno.loc[:,"Carcinogen_classification_annotation"] = dftoxianno.apply(Carcinogen_classification_annotation, axis = 1)
    dftoxianno["Hepatotoxicity_final_annotation"] = ""
    dftoxianno["Carcinogenicity_final_annotation"] = ""
    dftoxianno["Hepatotoxicity_literature review"] = ""
    dftoxianno["Carcinogenicity_literature review"] = ""
    dftoxianno.to_csv(r"in-silico toxicity and drug interaction prediction\integration\Toxi_infor_annota_sum.csv")