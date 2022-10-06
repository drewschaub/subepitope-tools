import annotation_tools
import configparser
import csv
import pandas as pd
import warnings
from abnumber import Chain
from pathlib import Path
from tqdm import tqdm

def peptide_label_df(df):
    """Generates a new column/feature with peptide description (e.g. heavy, light, fusion peptide, spike, gp41)

    Args:
        df: an existing dataframe that needs column names: assignment_name, pdbx_description, pdbx_seq_one_letter_code, 
    """    
    pdbx_description = df['pdbx_description'].lower()
    try:
        assignment_name = df['assignment_name'].lower()
    except:
        assignment_name = df['assignment_name']
    
    heavy_substrings = ['heavy', 'vh', 'vh-ch1']
    light_substrings = ['light', 'vl', 'kappa', 'lambda']
    gp160_substrings = ['gp160', 'glycoprotein 160']
    gp140_substrings = ['gp140', 'glycoprotein 140']
    gp120_substrings = ['gp120', 'glycoprotein 120']
    gp41_substrings = ['gp41', 'glycoprotein 41']
    fusion_peptide_substrings = ['fusion peptide']
    nanobody_substrings = ['nanobody', 'j3 vhh']
    tcell_surface_gp_cd4_substrings = ['cell surface glycoprotein cd4']
    adnectin_substrings = ['adnectin']
    M48U1_substrings = ['m48u1']
    
    if 'AVGIGAVF' in df['pdbx_seq_one_letter_code']:
        return 'fusion peptide'
    elif any(substring in pdbx_description for substring in nanobody_substrings):
        return 'nanobody'
    elif any(substring in pdbx_description for substring in heavy_substrings):
        return 'heavy'
    elif any(substring in pdbx_description for substring in light_substrings):
        return 'light'
    elif any(substring in pdbx_description for substring in gp160_substrings):
        return 'gp160'
    elif any(substring in pdbx_description for substring in gp140_substrings):
        return 'gp140'
    elif any(substring in pdbx_description for substring in gp120_substrings):
        return 'gp120'
    elif any(substring in pdbx_description for substring in gp41_substrings):
        return 'gp41'
    elif any(substring in pdbx_description for substring in fusion_peptide_substrings):
        return 'fusion peptide'
    elif any(substring in pdbx_description for substring in tcell_surface_gp_cd4_substrings):
        return 't-cell surface glycoprotein cd4'
    elif any(substring in pdbx_description for substring in adnectin_substrings):
        return 'adnectin'
    elif any(substring in pdbx_description for substring in M48U1_substrings):
        return 'm48u1'
    # Check in assignmen_name column
    if any(substring in assignment_name for substring in nanobody_substrings):
        return 'nanobody'
    elif any(substring in assignment_name for substring in heavy_substrings):
        return 'heavy'
    elif any(substring in assignment_name for substring in light_substrings):
        return 'light'
    elif any(substring in assignment_name for substring in gp160_substrings):
        return 'gp160'
    elif any(substring in assignment_name for substring in gp140_substrings):
        return 'gp140'
    elif any(substring in assignment_name for substring in gp120_substrings):
        return 'gp120'
    elif any(substring in assignment_name for substring in gp41_substrings):
        return 'gp41'
    elif any(substring in assignment_name for substring in fusion_peptide_substrings):
        return 'fusion peptide'
    elif any(substring in assignment_name for substring in tcell_surface_gp_cd4_substrings):
        return 't-cell surface glycoprotein cd4'
    elif any(substring in assignment_name for substring in adnectin_substrings):
        return 'adnectin'
    elif any(substring in assignment_name for substring in M48U1_substrings):
        return 'm48u1'
    else:
        return 'None Detected'

def main():
    # Load the config file
    config = configparser.ConfigParser()
    config.read(Path('..','..','config', 'config.ini'))

    # Get parameters from config file
    data_root = config['paths']['data_root']
    project_root = config['paths']['project_root']
    target = config['target']['antigen']
    api_url = config['paths']['rcsb_api_url']
    file_prefix = config['paths']['file_prefix']

    # Check if flag is set to supress warnings
    if config['misc']['ignore_warnings'] == 'yes':
        warnings.filterwarnings('ignore')

    # Check if flag is set to supress warnings
    if config['misc']['verbose'] == 'yes':
        verbose = True
    else:
        verbose = False

    # Set filetype paths and make directories if necessary 
    project_path = Path(project_root)
    data_path = Path(data_root, 'data', target)
    cif_path = Path(data_path, 'cif')
    cif_path.mkdir(parents=True, exist_ok=True)
    csv_path = Path(data_path, 'csv')
    csv_path.mkdir(exist_ok=True)
    fasta_path = Path(data_path, 'fasta')
    fasta_path.mkdir(exist_ok=True)
    json_path = Path(data_path, 'json')
    json_path.mkdir(exist_ok=True)
    pdb_path = Path(data_path, 'pdb')
    pdb_path.mkdir(exist_ok=True)
    pickle_path = Path(data_path, 'pickle')
    pickle_path.mkdir(exist_ok=True)
    txt_path = Path(data_path, 'txt')
    txt_path.mkdir(exist_ok=True)

    print('### Loading all CIF files')
    # All of the CIF files in our cif folder
    cif_paths = sorted(cif_path.glob('*.cif'))
    PDBid_list = annotation_tools.generate_stem_list(cif_paths)

    if verbose:
        print('{} unique PDB IDs'.format(len(PDBid_list)))

    # create an Empty DataFrame object
    df = pd.DataFrame(columns = ['pdbx_description',
                                'rcsb_macromolecular_names_combined_name',
                                'PDB',
                                'pdbx_strand_id',
                                'pdbx_strand_id_count', 
                                'polymer_entity_id',
                                'rcsb_sample_sequence_length',
                                'pdbx_seq_one_letter_code',
                                'pdbx_gene_src_scientific_name',
                                'provenance_code',
                                'reference_database_name',
                                'reference_database_accession',
                                'assignment_feature_id',
                                'assignment_name',
                                'assignment_provenance_source',
                                'assignment_type',
                                'URL',
                                'json_filename'
                                ])

    print('\n### Processing JSON files')
    if config['misc']['load_existing_json'] != 'yes':
        # fetch json files
        fetch_jsons(PDBid_list, json_path, target)

        # save the dataframe to a pickle and csv object for re-usability
        df.to_pickle(Path(pickle_path, '{}_polymer_entity.pkl'.format(target)))
        df.to_csv(Path(csv_path, '{}_polymer_entity.csv'.format(target)))
    else:
        df_working = pd.read_pickle(Path(pickle_path, '{}_polymer_entity.pkl'.format(target)))
        if verbose:
            print('{} records found in {}_polymer_entity.pkl'.format(len(df_working), target))
        df_working_PDB_list = df_working['PDB'].tolist()

    print('\n\tScientific names present:        ')
    print('=====================================')
    print(df_working['pdbx_gene_src_scientific_name'].value_counts())

    df_working['peptide_label'] = df_working.apply(peptide_label_df, axis = 1)

    print('\nSaving {}_polymer_entity_peptide-labels.csv'.format(target))
    df_working.to_csv(Path(csv_path, '{}_polymer_entity_peptide-labels.csv'.format(target)))

    print('\n\tPeptipe labels addded:        ')
    print('=====================================')
    print(df_working['peptide_label'].value_counts())

    print('\n### Reading antibody assignments from PDBid_antibodyname.csv')
    # User will need to create their own PDBid_antibodyname.csv
    reader = csv.DictReader(open(Path(csv_path, "PDBid_antibodyname.csv")))
    pdb_to_antibody_label_dict = {}
    unique_antibody_list = []

    for row in reader:
        key = row['PDBid']
        value = row['antibody_name']
        stored_value_list = []
        updated_value = ''
        
        if value not in unique_antibody_list:
            unique_antibody_list.append(value)
        
        if key in pdb_to_antibody_label_dict:
            stored_value_list = pdb_to_antibody_label_dict[key].split(';')
    
            if value not in stored_value_list:
                updated_value = '{};{}'.format(pdb_to_antibody_label_dict[key], value)
                pdb_to_antibody_label_dict[key] = updated_value
        else:
            pdb_to_antibody_label_dict[key] = value

    if verbose:
        print('{} unique antibody names detected'.format(len(unique_antibody_list)))
        print('{} unique entries in PDB:antibody_name dictionary due to some PDB structures having more than one antibody bound'.format(len(pdb_to_antibody_label_dict)))

    print('\n### Generating csv of unique PDBs with antibodies present as unique_PDB_antibodies.csv')
    # Generate a series that is unique PDBs, so may have multiple antibodies
    s = pd.Series(pdb_to_antibody_label_dict)
    s.index.name = 'PDB'
    unique_PDB_antibodies_df = s.to_frame('antibodies')
    unique_PDB_antibodies_df.to_csv(Path(csv_path, "unique_PDB_antibodies.csv"))

    print('\n### Generating csv of unique antibodies and which PDBs they are present as unique_antibody_PDBs.csv')
    # Generate a series that is unique antibodies, so may have multiple PDBs
    antibody_label_to_pdb_dict = {}
    for antibody in unique_antibody_list:

        reader = csv.DictReader(open(Path(csv_path, "PDBid_antibodyname.csv")))
        for row in reader:
            value = row['PDBid']
            key = row['antibody_name']
            stored_value_list = []
            updated_value = ''
        
            if key not in antibody_label_to_pdb_dict:
                antibody_label_to_pdb_dict[key] = value
            else:
                stored_value_list = antibody_label_to_pdb_dict[key].split(';')

                if value not in stored_value_list:
                    updated_value = '{};{}'.format(antibody_label_to_pdb_dict[key], value)
                    antibody_label_to_pdb_dict[key] = updated_value

    # Generate a series that is unique antibodies, so may have multiple PDBs
    s = pd.Series(antibody_label_to_pdb_dict)
    s.index.name = 'antibody'
    unique_antibody_PDBs_df = s.to_frame('PDBs')
    unique_antibody_PDB_list = unique_antibody_PDBs_df['PDBs'].tolist()
    PDB_count = []
    for PDBlist in unique_antibody_PDB_list:
        PDB_count.append(len(PDBlist.split(';')))
    unique_antibody_PDBs_df['PDB_count'] = PDB_count 
    unique_antibody_PDBs_df.to_csv(Path(csv_path, "unique_antibody_PDBs.csv"))

    if verbose:
        print('{} unique entries in antibody_name:PDB dictionary'.format(len(antibody_label_to_pdb_dict)))

    print('\n### Appending antibody label to main dataframe {}_polymer_entity'.format(target))
    reader = csv.DictReader(open(Path(csv_path, "unique_PDB_antibodies.csv")))
    pdb_to_antibody_label_dict = {}

    for row in reader:
        pdb_to_antibody_label_dict[row['PDB']] = row['antibodies']

    df_working['antibodies_in_PDB'] = df_working['PDB'].map(pdb_to_antibody_label_dict)

    new_antibody_labels = []
    for index, row in df_working.iterrows():
        if ';' in str(row['antibodies_in_PDB']):
            if row['peptide_label'] in ['heavy', 'light']:
                antibody_names = row['antibodies_in_PDB'].split(';')
                not_detected = True
                number_detected = 0
                for name in antibody_names:
                    if name in row['pdbx_description']:
                        number_detected = number_detected + 1
                        new_antibody_label = name
                        not_detected = False
                        if number_detected > 1:
                            print('WARNING!!, More than one antibody from: {} matches whiel parsing {} in PDB: {}'.format(antibody_names, row['pdbx_description'], row['PDB']))
                            new_antibody_label = 'more-than-one-match'
                if not_detected:
                    if (row['PDB'] == "7PC2") and (row['pdbx_strand_id'] == "K,G,I") and (row['peptide_label'] == "heavy"):
                        name = "3BNC117"
                    else:
                        print('WARNING!! No antibodies from: {} was not detected while parsing {} in PDB: {}'.format(antibody_names, row['pdbx_description'], row['PDB']))
                        new_antibody_label = 'no-match'
                        not_detected = True
                    
                new_antibody_labels.append(new_antibody_label)
            else:
                new_antibody_labels.append(row['antibodies_in_PDB'])
        else:
            new_antibody_labels.append(row['antibodies_in_PDB'])

    df_working['antibody_label'] = new_antibody_labels

    if config['misc']['use_class_definition_csv'] == 'yes':
        df_class_defs = pd.read_csv(Path(csv_path, 'class_definitions.csv'))
        df_working['class'] = df_working['antibody_label']
        df_working['category'] = df_working['antibody_label']
        antibody_class_dict = dict(zip(df_class_defs['antibody'], df_class_defs['class']))
        antibody_category_dict = dict(zip(df_class_defs['antibody'], df_class_defs['category']))

        df_working['class'] = df_working['class'].map(antibody_class_dict)
        df_working['category'] = df_working['category'].map(antibody_category_dict)
        print(df_working)

    # save the dataframe to a pickle and csv object for re-usability
    df.to_csv(Path(csv_path, '{}_polymer_entity_labeled.csv'.format(target)))

    if config['misc']['add_ANARCI_data'] == 'yes':
        sequence_list = df_working['pdbx_seq_one_letter_code'].tolist()

        chain_types = []
        cdr1_seqs = []
        cdr2_seqs = []
        cdr3_seqs = []
        cdr_defs = []
        fr1_seqs = []
        fr2_seqs = []
        fr3_seqs = []
        fr4_seqs = []
        j_genes = []
        v_genes = []
        names = []
        schemes = []
        species = []
        tails = []
        is_heavy = []
        is_light = []
        is_kappa = []
        is_lambda = []

        for seq in tqdm(sequence_list):
            
            try:
                try:
                    chain = Chain(seq, scheme='kabat', assign_germline=True)
                except:
                    chain = Chain(seq, scheme='kabat')
                    print('Failed to assign germline')

                try:
                    chain_types.append(chain.chain_type)
                except:
                    chain_types.append(None)
                    
                try:
                    cdr1_seqs.append(chain.cdr1_seq)
                except:
                    cdr1_seqs.append(None)
                
                try:
                    cdr2_seqs.append(chain.cdr2_seq)
                except:
                    cdr2_seqs.append(None)
                    
                try:
                    cdr3_seqs.append(chain.cdr3_seq)
                except:
                    cdr3_seqs.append(None)
                    
                try:
                    cdr_defs.append(chain.cdr_definition)
                except:
                    cdr_defs.append(None)
                    
                try:
                    fr1_seqs.append(chain.fr1_seq)
                except:
                    fr1_seqs.append(None)
                
                try:
                    fr2_seqs.append(chain.fr2_seq)
                except:
                    fr2_seqs.append(None)
                    
                try:
                    fr3_seqs.append(chain.fr3_seq)
                except:
                    fr3_seqs.append(None)

                try:
                    fr4_seqs.append(chain.fr4_seq)
                except:
                    fr4_seqs.append(None)
                    
                try:
                    j_genes.append(chain.j_gene)
                except:
                    j_genes.append(None)
                    
                try:
                    v_genes.append(chain.v_gene)
                except:
                    v_genes.append(None)
                    
                try:
                    names.append(chain.name)
                except:
                    names.append(None)
                    
                try:
                    schemes.append(chain.scheme)
                except:
                    schemes.append(None)
                    
                try:
                    species.append(chain.species)
                except:
                    species.append(None)
                    
                try:
                    tails.append(chain.tail)
                except:
                    tails.append(None)
                    
                try:
                    is_heavy.append(chain.is_heavy_chain())
                except:
                    is_heavy.append(False)
                    
                try:
                    is_light.append(chain.is_light_chain())
                except:
                    is_light.append(False)
                    
                try:
                    is_kappa.append(chain.is_kappa_chain())
                except:
                    is_kappa.append(False)
                    
                try:
                    is_lambda.append(chain.is_lambda_chain())
                except:
                    is_lambda.append(False)
                    
            except:
                chain_types.append(None)
                is_heavy.append(False)
                is_light.append(False)
                is_kappa.append(False)
                is_lambda.append(False)  
                cdr1_seqs.append(None)
                cdr2_seqs.append(None)
                cdr3_seqs.append(None)
                cdr_defs.append(None)
                fr1_seqs.append(None)
                fr2_seqs.append(None)
                fr3_seqs.append(None)
                fr4_seqs.append(None)
                j_genes.append(None)
                v_genes.append(None)
                names.append(None)
                schemes.append(None)
                species.append(None)
                tails.append(None)
                
        df_working['chain_type_ANARCI'] = chain_types
        df_working['is_heavy_ANARCI'] = is_heavy
        df_working['is_light_ANARCI'] = is_light
        df_working['is_kappa_ANARCI'] = is_kappa
        df_working['is_lambda_ANARCI'] = is_lambda
        df_working['cdr1_seq_ANARCI'] = cdr1_seqs
        df_working['cdr2_seq_ANARCI'] = cdr2_seqs
        df_working['cdr3_seq_ANARCI'] = cdr3_seqs
        df_working['cdr_def_ANARCI'] = cdr_defs
        df_working['fr1_seq_ANARCI'] = fr1_seqs
        df_working['fr2_seq_ANARCI'] = fr2_seqs
        df_working['fr3_seq_ANARCI'] = fr3_seqs
        df_working['fr4_seq_ANARCI'] = fr4_seqs
        df_working['j_gene_ANARCI'] = j_genes
        df_working['v_gene_ANARCI'] = v_genes
        df_working['name_ANARCI'] = names
        df_working['scheme_ANARCI'] = schemes
        df_working['species_ANARCI'] = species
        df_working['tail_ANARCI'] = tails

        df_working.to_pickle(Path(pickle_path, '{}_polymer_entity_ANARCI.pkl'.format(target)))
        columns = ['PDB',
                'pdbx_strand_id',
                'pdbx_description',
                'rcsb_sample_sequence_length',
                'antibodies_in_PDB',
                'antibody_label',
                'category',
                'class',
                'peptide_label',
                'pdbx_gene_src_scientific_name',
                'species_ANARCI',
                'chain_type_ANARCI',
                'is_heavy_ANARCI',
                'is_light_ANARCI',
                'is_kappa_ANARCI',
                'is_lambda_ANARCI',
                'scheme_ANARCI',
                'cdr_def_ANARCI',
                'fr1_seq_ANARCI',
                'cdr1_seq_ANARCI',
                'fr2_seq_ANARCI',
                'cdr2_seq_ANARCI',
                'fr3_seq_ANARCI',
                'cdr3_seq_ANARCI',
                'fr4_seq_ANARCI',
                'v_gene_ANARCI',
                'j_gene_ANARCI', 
                'tail_ANARCI',
                'name_ANARCI',
                'reference_database_name',
                'reference_database_accession',
                'assignment_feature_id',
                'assignment_name',
                'assignment_provenance_source',
                'assignment_type',
                'provenance_code',
                'URL',
                'json_filename',
                'pdbx_strand_id_count',
                'polymer_entity_id',
                'rcsb_macromolecular_names_combined_name', 
                'pdbx_seq_one_letter_code']

        df_working = df_working[columns]
        df_working.to_csv(Path(csv_path, '{}_polymer_entity_ANARCI.csv'.format(target)))       

        antigen_csv_filename = '{}_antigen_{}_PDBscrape.csv'.format(file_prefix,target)
        h_chain_csv_filename = '{}_h_chain_{}_PDBscrape.csv'.format(file_prefix,target)
        l_chain_csv_filename = '{}_l_chain_{}_PDBscrape.csv'.format(file_prefix,target)
        antibody_csv_filename = '{}_antibody_{}_PDBscrape.csv'.format(file_prefix,target)

        print('\n### Verifying heavy and light chain assignments using ANARCI')
        print('\n### Saving summary files for antigen, heavy chain and light chain as:')
        print('      - {}'.format(antigen_csv_filename))
        print('      - {}'.format(h_chain_csv_filename))
        print('      - {}'.format(l_chain_csv_filename))
        print('      - {}'.format(antibody_csv_filename))

        df_antigen = df_working[df_working['peptide_label'] == 'fusion peptide']
        df_h_chain = df_working[df_working['is_heavy_ANARCI'] == True]
        df_l_chain = df_working[df_working['is_light_ANARCI'] == True]
        df_antibodies = df_h_chain.append(df_l_chain)

        df_antigen.to_csv(Path(csv_path, antigen_csv_filename))
        df_h_chain.to_csv(Path(csv_path, h_chain_csv_filename))
        df_l_chain.to_csv(Path(csv_path, l_chain_csv_filename))
        df_antibodies.to_csv(Path(csv_path, antibody_csv_filename))

        print('\n### Heavy and light chain summary')
        print('        {} unique light chains'.format(len(df_l_chain['pdbx_seq_one_letter_code'].unique().tolist())))
        print('        {} unique heavy chains'.format(len(df_h_chain['pdbx_seq_one_letter_code'].unique().tolist())))

    final_structure_summary_csvfile = '{}_{}_structure_summary.csv'.format(file_prefix, target)
    df_working = df_working.sort_values(by=['peptide_label', 'antibody_label', 'PDB']).reset_index(drop=True)
    df_working.to_csv(Path(csv_path, final_structure_summary_csvfile))

    # create an dataframe that will be used to generate file stems with antibody label, peptide label, PDB id, chain id
    df_pdb_assignments = pd.DataFrame(columns = ['antibody_label', 'peptide_label', 'pdb', 'chain', 'file_prefix'])

    for index, row in df_working.iterrows():
        antibody_label = row['antibody_label']
        peptide_label = row['peptide_label']
        pdb = row['PDB']
        
        temp_list = row['pdbx_strand_id'].split(',')
        chain_list = []
        for chain in temp_list:
            chain_list.append(chain[0])
            if len(chain) > 1:
                print('WARNING: Chain {} is more than one character in PDB: {}. Not allowed in PDB, but allowed in CIF. PDB may not exist on RCSB PDB'.format(chain, pdb))
        
        for chain in chain_list:
            file_prefix = "{}_{}_{}_{}".format(antibody_label, peptide_label, pdb, chain)
            polymer_entity_dict = {'antibody_label' : antibody_label, 'peptide_label' : peptide_label, 'pdb' : pdb,'chain' : chain, 'file_prefix' : file_prefix}
            df_pdb_assignments = df_pdb_assignments.append(polymer_entity_dict, ignore_index=True)

    pdb_assignment_filename = "PDB_chain_peptide-label_antibody-label_file-prefix.csv"
    print('\n### Generating csv file that will contain file naming structure')
    print('      PDB id, chain id, peptide label, antibody label, file prefix')
    print('      {}'.format(pdb_assignment_filename))
    df_pdb_assignments.to_csv(Path(csv_path, pdb_assignment_filename))

    print('\n### Double Check Structure Summary File')
    print('')
    print('        {}'.format(final_structure_summary_csvfile))


if __name__ == "__main__":
    main()