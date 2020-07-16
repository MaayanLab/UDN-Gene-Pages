#!/usr/bin/env python
# coding: utf-8

# In[12]:


#! /usr/bin/python

import json
import requests
import pandas as pd
from collections import defaultdict


#read entire database into dataframe

#create list for final table
genes = []
matches = []
enrichr = []
enrichr_hit = []
#My data 
my_data = 'https://docs.google.com/spreadsheets/d/1N038nU0Xc1_Ob2cA97HzjbgqE4KCdKbGb7VahJtJzsY/edit#gid=0'


#Turn CSV export URL for my data
csv_export_url = my_data.replace('/edit#gid=', '/export?format=csv&gid=')

#df is entire UDN database
df = pd.read_csv(csv_export_url)
UDN_dict = df.set_index('Gene').groupby(level=0).agg(list).to_dict()['HGO_Disease_ID']

#print(UDN_dict['PPFIA3'])

#returns all the keys
def getList(dict): 
    return dict.keys()

# function to return key for any value 
def get_key(val): 
    for key, value in my_dict.items(): 
         if val == value: 
            return key 
  
    return "key doesn't exist"


#iterate through each key

print("ENRICHR PREDICTION METHOD MACTHES: ")
for key, value in UDN_dict.items():
    #keep track of phenotypes from my table and add to a new list
    phenotype = value
    #keep track of gene
    gene = key
    #print(gene)
    #phenotype.append(value)
    #print(value)
    
    # ENRICHR DATA
    URL = 'https://amp.pharm.mssm.edu/geneshot/api/predict/gene/'
    BACK_URL = '/library/Human_Phenotype_Ontology/similarity/enrichr/offset/0/limit/200'
    ENRICHR_URL = URL + gene + BACK_URL 
    
    ENRICHR_response = requests.get(
    ENRICHR_URL
 )
    if not ENRICHR_response.ok:
        raise Exception('Error during query')
    
    ENRICHR_data = json.loads(ENRICHR_response.text)
    
    ENRICHR_results = ENRICHR_data.get('results')
    
    #if result is empty, this was not found on geneshot so break 
    if ENRICHR_results is None:
        continue
    

    
    #store phenotypes from geneshot data
    a = []
    rank = 1
    
    for dict in ENRICHR_results:
        rank+=1
        for key,value in dict.items():
            a.append(dict['property'])
            break
    #print(a)
    
    #find matches
    x = set(a) & set(phenotype)
    
    
    #create a list of matches with ranks
    b = []
    w = 0;
    for val in x:
        w = a.index(val);
        w+=1;
        z= 'Rank:' + str(w)
        b.append([val, w])  
    
    #sort the list by ranks
    b.sort(key=lambda x:x[1])
    
    #final table
    genes.append(gene)
    matches.append(b)
    enrichr.append('Enrichr')
    enrichr_hit.append(len(b))
    


    print('Gene:' + gene)
    
    if not b:
          print("No matches found")
    for val in b:
        print(val)

    print()

#table for matches lists
autorif_gene = []
autorif_match = []
autorif = []
autorif_hit = []

#AUTORIF CODE
print("--------------------------------------------------------------------------")
print("AUTORIF PREDICTION METHOD MACTHES: ")
for key, value in UDN_dict.items():
    #keep track of phenotypes from my table and add to a new list
    phenotype = value
    #keep track of gene
    gene = key
    #print(gene)
    #phenotype.append(value)
    #print(value)
    
    # ENRICHR DATA
    URL = 'https://amp.pharm.mssm.edu/geneshot/api/predict/gene/'
    BACK_URL = '/library/Human_Phenotype_Ontology/similarity/autorif/offset/0/limit/200'
    ENRICHR_URL = URL + gene + BACK_URL 
    
    ENRICHR_response = requests.get(
    ENRICHR_URL
 )
    if not ENRICHR_response.ok:
        raise Exception('Error during query')
    
    ENRICHR_data = json.loads(ENRICHR_response.text)
    
    ENRICHR_results = ENRICHR_data.get('results')
    
    #if result is empty, this was not found on geneshot so break 
    if ENRICHR_results is None:
        continue
    

    
    #store phenotypes from geneshot data
    a = []
    rank = 1
    
    for dict in ENRICHR_results:
        rank+=1
        for key,value in dict.items():
            a.append(dict['property'])
            break
    #print(a)
    
    #find matches
    x = set(a) & set(phenotype)
    
    
    #create a list of matches with ranks
    b = []
    w = 0;
    for val in x:
        w = a.index(val);
        w+=1;
        z= 'Rank:' + str(w)
        b.append([val, w])  
    
    #sort the list by ranks
    b.sort(key=lambda x:x[1])

    autorif_gene.append(gene)
    autorif_match.append(b)
    autorif.append('Autorif')
    autorif_hit.append(len(b))
    print('Gene:' + gene)
    #final table
   
    
    if not b:
          print("No matches found")
    for val in b:
        print(val)

    print()

#create list for final table generif
generif_genes = []
generif_matches = []
generif = []
generif_hit = []


#GENERIF CODE
print("--------------------------------------------------------------------------")
print("GENERIF PREDICTION METHOD MACTHES: ")
for key, value in UDN_dict.items():
    #keep track of phenotypes from my table and add to a new list
    phenotype = value
    #keep track of gene
    gene = key
    #print(gene)
    #phenotype.append(value)
    #print(value)
    
    # ENRICHR DATA
    URL = 'https://amp.pharm.mssm.edu/geneshot/api/predict/gene/'
    BACK_URL = '/library/Human_Phenotype_Ontology/similarity/generif/offset/0/limit/200'
    ENRICHR_URL = URL + gene + BACK_URL 
    
    ENRICHR_response = requests.get(
    ENRICHR_URL
 )
    if not ENRICHR_response.ok:
        raise Exception('Error during query')
    
    ENRICHR_data = json.loads(ENRICHR_response.text)
    
    ENRICHR_results = ENRICHR_data.get('results')
    
    #if result is empty, this was not found on geneshot so break 
    if ENRICHR_results is None:
        continue
    

    
    #store phenotypes from geneshot data
    a = []
    rank = 1
    
    for dict in ENRICHR_results:
        rank+=1
        for key,value in dict.items():
            a.append(dict['property'])
            break
    #print(a)
    
    #find matches
    x = set(a) & set(phenotype)
    
    
    #create a list of matches with ranks
    b = []
    w = 0;
    for val in x:
        w = a.index(val);
        w+=1;
        z= 'Rank:' + str(w)
        b.append([val, w])  
    
    generif_genes.append(gene)
    generif_matches.append(b)
    generif.append('generif')
    generif_hit.append(len(b))
    
    #sort the list by ranks
    b.sort(key=lambda x:x[1])

    print('Gene:' + gene)
    
    if not b:
          print("No matches found")
    for val in b:
        print(val)

    print()
    
#create list for final table tagger
tagger_genes = []
tagger_matches = []
tagger = []
tagger_hit = []

#TAGGER CODE
print("--------------------------------------------------------------------------")
print("TAGGER PREDICTION METHOD MACTHES: ")
for key, value in UDN_dict.items():
    #keep track of phenotypes from my table and add to a new list
    phenotype = value
    #keep track of gene
    gene = key
    #print(gene)
    #phenotype.append(value)
    #print(value)
    
    # ENRICHR DATA
    URL = 'https://amp.pharm.mssm.edu/geneshot/api/predict/gene/'
    BACK_URL = '/library/Human_Phenotype_Ontology/similarity/tagger/offset/0/limit/200'
    ENRICHR_URL = URL + gene + BACK_URL 
    
    ENRICHR_response = requests.get(
    ENRICHR_URL
 )
    if not ENRICHR_response.ok:
        raise Exception('Error during query')
    
    ENRICHR_data = json.loads(ENRICHR_response.text)
    
    ENRICHR_results = ENRICHR_data.get('results')
    
    #if result is empty, this was not found on geneshot so break 
    if ENRICHR_results is None:
        continue
    

    
    #store phenotypes from geneshot data
    a = []
    rank = 1
    
    for dict in ENRICHR_results:
        rank+=1
        for key,value in dict.items():
            a.append(dict['property'])
            break
    #print(a)
    
    #find matches
    x = set(a) & set(phenotype)
    
    
    #create a list of matches with ranks
    b = []
    w = 0;
    for val in x:
        w = a.index(val);
        w+=1;
        z= 'Rank:' + str(w)
        b.append([val, w])  
    
    tagger_genes.append(gene)
    tagger_matches.append(b)
    tagger.append('tagger')
    tagger_hit.append(len(b))
    
    #sort the list by ranks
    b.sort(key=lambda x:x[1])

    print('Gene:' + gene)
    
    if not b:
          print("No matches found")
    for val in b:
        print(val)

    print()

#create list for final table tagger
archs4_genes = []
archs4_matches = []
archs4= []
archs4_hit = []
column_names = ['Gene','Matches from Archs4']

dff = pd.DataFrame(columns = column_names)

#ARCHS4 CODE
print("--------------------------------------------------------------------------")
print("ARCHS4 PREDICTION METHOD MACTHES: ")
for key, value in UDN_dict.items():
    #keep track of phenotypes from my table and add to a new list
    phenotype = value
    #keep track of gene
    gene = key
    #print(gene)
    #phenotype.append(value)
    #print(value)
    
    # ENRICHR DATA
    URL = 'https://amp.pharm.mssm.edu/geneshot/api/predict/gene/'
    BACK_URL = '/library/Human_Phenotype_Ontology/similarity/coexpression/offset/0/limit/200'
    ENRICHR_URL = URL + gene + BACK_URL
    
    ENRICHR_response = requests.get(
    ENRICHR_URL
 )
    if not ENRICHR_response.ok:
        raise Exception('Error during query')
    
    ENRICHR_data = json.loads(ENRICHR_response.text)
    
    ENRICHR_results = ENRICHR_data.get('results')
    
    #if result is empty, this was not found on geneshot so break 
    if ENRICHR_results is None:
        continue
    

    
    #store phenotypes from geneshot data
    a = []
    rank = 1
    
    for dict in ENRICHR_results:
        rank+=1
        for key,value in dict.items():
            a.append(dict['property'])
            break
    #print(a)
    
    #find matches
    x = set(a) & set(phenotype)
    
    
    #create a list of matches with ranks
    b = []
    w = 0;
    for val in x:
        w = a.index(val);
        w+=1;
        z= 'Rank:' + str(w)
        b.append([val, w])  
    
   
        
    length = len(b)
    print(len(b))
    dff = pd.DataFrame({'Gene': genes, 'Matches from Archs4': length })

    #sort the list by ranks
    b.sort(key=lambda x:x[1])
    
    archs4_genes.append(gene)
    archs4_matches.append(b)
    archs4.append('archs4')
    archs4_hit.append(len(b))
    

    print('Gene:' + gene)
    
    if not b:
          print("No matches found")
    for val in b:
        print(val)

    print()

final = pd.DataFrame({'Matches': matches, 'Gene': genes, 'Type': enrichr, '#Enrichr Hits': enrichr_hit})

autorif_final = pd.DataFrame({'Matches': autorif_match, 'Gene': autorif_gene, 'Type': autorif, '#Autorif Hits': autorif_hit})

generif_final = pd.DataFrame({'Matches': generif_matches, 'Gene': generif_genes, 'Type': generif,'#Autorif Hits': generif_hit})
tagger_final = pd.DataFrame({'Matches': tagger_matches, 'Gene': tagger_genes, 'Type': tagger, '#Tagger Hits': tagger_hit})
archs4_final = pd.DataFrame({'Matches': archs4_matches, 'Gene': archs4_genes, 'Type': archs4,'#ARCHS4 Hits': archs4_hit})

#merge different datasets

new = pd.merge(final, autorif_final, on="Gene", how = "right")
news = pd.merge(new, generif_final, on="Gene",  how = "right")
tagger= pd.merge(news, tagger_final, on="Gene",  how = "right")
archs4 = pd.merge(tagger, archs4_final, on="Gene",  how = "right")




# In[13]:


print(archs4)


# In[15]:


archs4.to_csv (r'C:\Users\samanthabuuiyan\Desktop\HGO_TABLE.csv', index = False, header=True)


# In[ ]:




