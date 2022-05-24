import pandas as pd
import numpy as np
import sys
from pronto import Ontology
   

def create_agreement_cols(df_to_work_filled_merge:pd.DataFrame) -> pd.DataFrame:
    """Receives a df and save a new one as
    a csv file including the Diff columns
    (agree/disagree with ontology description)"""

    df_to_work_filled_merge['Diff_health'] = np.where(df_to_work_filled_merge['donor_health_status'] == df_to_work_filled_merge['Descrip_onto_health_status'] , 'agree', 'disagree')
    df_to_work_filled_merge['Diff_health_merge'] = np.where(df_to_work_filled_merge['donor_health_status_merge'] == df_to_work_filled_merge['Descrip_onto_health_status'] , 'agree', 'disagree')
    df_to_work_filled_merge['Diff_disease'] = np.where(df_to_work_filled_merge['disease'] == df_to_work_filled_merge['Descrip_onto_disease'] , 'agree', 'disagree')
    
   #reorder columns
    cols = ['EpiRR', 'EpiRR_status', 'age', 'biomaterial_type', 'cell_type',
       'disease_ontology_term', 'donor_age_unit',
       'donor_id', 'donor_life_stage', 'health_state', 'line', 'markers',
       'project', 'sample_ontology_term', 'sex', 'taxon_id', 'tissue_type',
       'donor_health_status_merge', 'donor_health_status', 'Descrip_onto_health_status', 
       'Diff_health', 'Diff_health_merge', 'donor_health_status_ontology_uri', 
       'donor_health_status_ontology_curie','disease', 'Descrip_onto_disease', 'Diff_disease',
       'disease_ontology_uri', 'disease_ontology_curie', 'Descrip_all_terms']
    
    df_to_work_filled_merge = df_to_work_filled_merge[cols]
    df_to_work_filled_merge.to_csv('IHEC_diff_onto_disease_health_ncitowl_2.csv', index=False) #saving file to explore


def fill_health_disease(epirr_amed_113,df_to_work:pd.DataFrame) -> pd.DataFrame:
    """Receives a list of EpiRR and a df with 
    health and health merge column. Returns a df
    """

    df_to_work_filled_merge = fill_amed_crest_merge(epirr_amed_113, df_to_work)

    #Fill na disease cells with descrp_disease
    df_to_work_filled_merge['disease'] = np.where((pd.isna(df_to_work_filled_merge['disease']) & \
        pd.notna(df_to_work_filled_merge['Descrip_onto_disease'])), 
        df_to_work_filled_merge['Descrip_onto_disease'], 
        df_to_work_filled_merge['disease'])

    #Fill na heath status cells with descrp_health
    df_to_work_filled_merge['donor_health_status'] = np.where((pd.isna(df_to_work_filled_merge['donor_health_status']) & \
        pd.notna(df_to_work_filled_merge['Descrip_onto_health_status'])), 
        df_to_work_filled_merge['Descrip_onto_health_status'], 
        df_to_work_filled_merge['donor_health_status'])

    #Fill na health status cells with descrp_all when disease is not na
    df_to_work_filled_merge['donor_health_status'] = np.where((pd.isna(df_to_work_filled_merge['donor_health_status']) & \
        pd.notna(df_to_work_filled_merge['disease'])), 
        df_to_work_filled_merge['Descrip_all_terms'], 
        df_to_work_filled_merge['donor_health_status'])
    
    return df_to_work_filled_merge


def fill_amed_crest_merge(epirr_amed_113, df_to_work:pd.DataFrame) -> pd.DataFrame:
    """Receives a list of EpiRR and a df with 
    health and health merge column. Returns a df
    """

    #whole dict 
    dict_merge = dict(zip(df_to_work['EpiRR'], df_to_work['donor_health_status_merge']))

    #list to filter dict (113 - AMED-CREST)
    list_to_dict = [line.strip() for line in epirr_amed_113]

    #filtering dict merge to map health column with new info
    dict_to_map = {my_key:dict_merge[my_key] for my_key in list_to_dict}

    #map dict
    df_to_work['donor_health_status'] = df_to_work['EpiRR'].map(dict_to_map)

    return df_to_work
    
 
def map_term_ncit(df_to_work:pd.DataFrame, ncit_obo:dict, ncit_dat:dict) -> pd.DataFrame:
    """Receives a df with the desired columns
    and two dictionaries (from .obo and .dat
    files). Return a df including the description
    NCI terms columns."""
    
    dict_terms_nci = create_ncit_obo_dict(ncit_obo)
    dict_dat = create_dict_dat(ncit_dat)

    #list ontologies in disease_ontology_term merged
    dis_ont = df_to_work['disease_ontology_term'].str.split(':').str[-1].tolist() #problem merged ::

    #list ontologies in disease_ontology_uri
    dis_ont_uri = df_to_work['disease_ontology_uri'].str.split('code=').str[-1].str.split('&').str[0].to_list()
    df_to_work['disease_ontology_uri'] = dis_ont_uri #reassigning values
    
    #list ontologies in donor_health_status_uri 
    dhealth_ont_uri = df_to_work['donor_health_status_ontology_uri'].str.split('code=').str[-1].str.split('&').str[0].tolist()
    df_to_work['donor_health_status_ontology_uri'] = dhealth_ont_uri

    #list to be columns
    health_ont_desc = []
    disease_ont_desc = []
    general_ont_desc = [] #to create as well

    health_ont_desc = create_description_col(dis_ont, dhealth_ont_uri, dict_terms_nci, dict_dat, health_ont_desc)
    disease_ont_desc = create_description_col(dis_ont, dis_ont_uri, dict_terms_nci, dict_dat, disease_ont_desc)

    #list with all descriptions available for all terms in disease_ontology_term
    for ele in dis_ont:
        if pd.notna(ele):

            if ele not in dict_terms_nci.keys():
                if ele == 'C0277545':
                    general_ont_desc.append('Disease type AND/OR category unknown')
                else:
                    term = dict_dat.get(ele) #getting value (ncit term) for CUI
                    general_ont_desc.append(dict_terms_nci.get(term, "No description available"))
            else:
                general_ont_desc.append(dict_terms_nci.get(ele, "No description available"))
        else:
            general_ont_desc.append(np.nan)

    #creating columns
    df_to_work['Descrip_onto_health_status'] = health_ont_desc
    df_to_work['Descrip_onto_disease'] = disease_ont_desc
    df_to_work['Descrip_all_terms'] = general_ont_desc

    return df_to_work


def create_description_col(list_ont:list, list_desired_col:list, dict_terms_nci:dict, dict_dat:dict, list_to_col:list) -> list:
    """Inputs: list ont, list_desired_col, dicts
       Output: list_to_col. The two input lists have
       the terms related with disease_ont_term and health
       or disease uri colums. The dicts were generated from
       .obo and .dat files. The output list constains the 
       description of each term"""

    for d_ont,d_uri in zip(list_ont, list_desired_col):
        if pd.notna(d_ont) and pd.notna(d_uri): 

            if d_uri not in dict_terms_nci.keys():
                if d_uri == 'C0277545':
                    list_to_col.append('Disease type AND/OR category unknown')
                else:
                    term = dict_dat.get(d_uri) #getting value (ncit term) for CUI in .dat dict
                    list_to_col.append(dict_terms_nci.get(term, "No description available"))
            else:
                list_to_col.append(dict_terms_nci.get(d_uri, "No description available")) #not working to replace None  
        else:
            list_to_col.append(np.nan)

    return list_to_col


def merge_col_ori_version(df_raw:pd.DataFrame, df_v7:pd.DataFrame) -> pd.DataFrame: #ok
    """Receives two dfs and returns a
    merged df (e.g IHEC metadata raw 
    and v7)"""

    df_ori_col = get_desired_cols(df_raw) #v01 - correct
    df_v7.rename({'donor_health_status': 'donor_health_status_merge'}, axis=1, inplace=True) #v07
    df_to_work = df_v7.merge(df_ori_col, on='EpiRR' ,how='left') #merging ori columns

    return df_to_work


def get_desired_cols(df_raw:pd.DataFrame) -> pd.DataFrame: #original version
    """Receives a df and returns
    a df with specific
    columns"""

    return df_raw[['EpiRR','donor_health_status','disease','donor_health_status_ontology_curie',
    'donor_health_status_ontology_uri','disease_ontology_curie',
    'disease_ontology_uri']]


def create_dict_dat(ncit_dat) -> dict:
    """Receives a .dat files and
    returns a dict (Keys=CUI term;
    values=NCIT term to map on
    obo dict)"""

    #keys CUI; val NCIT term
    return {line.strip().split('|')[1] : line.strip().split('|')[0]  for line in ncit_dat}


def create_ncit_obo_dict(ncit_obo) -> dict:
    """Receives a .obo file and
    returns a dict (Keys=terms;
    values=description_term)"""

    return {term.id.split('#')[-1] : term.name for term in ncit_obo.terms()}


def generate_obo_file(owl_file):
    """Receives a OWL file and 
    returns a .obo file"""

    with open("ncit_2204.obo", "wb") as f:
        owl_file.dump(f, format="obo")


def main():

    print('Starting script....')
    # ncit_owl = Ontology("ontology/Thesaurus.owl") #just run one time
    # generate_obo_file(nci_ont) #just run one time
    ncit_obo = Ontology("ontology/ncit_2204.obo")
    ncit_dat = open(sys.argv[1], 'r') #file .dat downloaded from ftp server (NCI)
    df_raw = pd.read_csv(sys.argv[2]) #version 1 (EpiRR updated - no merged columns)
    df_v7 = pd.read_csv(sys.argv[3]) #version 7 (last version)
    epirr_amed_113 = open(sys.argv[4],'r')
    df_to_work = merge_col_ori_version(df_raw,df_v7)
    df_to_work_1 = map_term_ncit(df_to_work, ncit_obo, ncit_dat)
    df_to_work_filled_merge = fill_health_disease(epirr_amed_113,df_to_work_1)

    print('Saving df...')
    create_agreement_cols(df_to_work_filled_merge)
    print('Finished!')


if __name__ == "__main__":

    main()