import annotation_tools
import configparser
import pandas as pd
import warnings
from abnumber import Chain
from Bio.PDB import *
from pathlib import Path

def runANARCI(infile, scheme = "kabat", species="human"):
    proc = subprocess.Popen(["ANARCI", "--sequence", infile, "--restrict", "ig", "--scheme", scheme, "--outfile", "out", "--csv", "--use_species", species], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    print("ANARCI", "--sequence", infile, "--restrict", "ig", "--scheme", scheme, "--outfile", "out", "--csv", "--use_species", species)
    outs, errs = proc.communicate()

def main():
    oneletter_dict = {
        'ASP':'D',
        'GLU':'E',
        'ASN':'N',
        'GLN':'Q',
        'ARG':'R',
        'LYS':'K',
        'PRO':'P',
        'GLY':'G',
        'CYS':'C',
        'THR':'T',
        'SER':'S',
        'MET':'M',
        'TRP':'W',
        'PHE':'F',
        'TYR':'Y',
        'HIS':'H',
        'ALA':'A',
        'VAL':'V',
        'LEU':'L',
        'ILE':'I'
    }

    residues_ignore_list = ['1PE', '44E', 'ACT', 'BR', 'BMA', 'BTB', 'CA', ' CA', 'CAC', 'CIT', 'CL', 'DMS', 'EDO', 'FMT', 
                            'GOL', 'HOH', 'IMD', 'IOD', 'K', 'MAN', 'MG', 'MLI', 'NA', 'NAG', 'NI', 'NO3', 'PCA', 'PEG', 
                            'PG0', 'PG4', 'PG5', 'PG6', 'PGE', 'PO4', 'SCN', 'SO4', 'TAR', 'TRS', 'UNK']


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

    print('\n### Loading all PDB files')
    # All of the CIF files in our cif folder
    pdb_paths = sorted(pdb_path.glob('*.pdb'))    

    PDBid_list = annotation_tools.generate_stem_list(cif_paths)

    if verbose:
        print('{} unique PDB IDs'.format(len(PDBid_list)))

    # Make dictionary of PDB (entry ID) to antibody
    unique_PDB_antibodies_file = "unique_PDB_antibodies.csv"
    unique_PDB_antibodies_path = Path(csv_path, unique_PDB_antibodies_file)
    df_unique_PDB_antibodies = pd.read_csv(unique_PDB_antibodies_path)

    keys = df_unique_PDB_antibodies['PDB'].tolist()
    values = df_unique_PDB_antibodies['antibodies'].tolist()
    PDB_antibody_dict = dict(zip(keys, values))

    # Make dataframe with chain information
    structure_summary_file = '{}_{}_structure_summary.csv'.format(file_prefix, target)
    structure_summary_path = Path(csv_path, structure_summary_file)
    df_structure_summary = pd.read_csv(structure_summary_path, index_col=0)
    print(df_structure_summary.head())

    exceptions_list = []
    count = 0
    for PDB_id in PDBid_list:
        print(PDB_id, count / len(PDBid_list))
        count = count + 1
        structure_name = "name"
        PDB_file = "{}.pdb".format(PDB_id)
        file_path = Path(pdb_path, PDB_file)

        # Genearte the H and L chain fasta files from the raw PDB file
        parser = PDBParser()
        ppb = PPBuilder()
        structure = parser.get_structure(PDB_id, file_path)

        # Trim dataframe to only current PDB id
        current_PDB_df = df_structure_summary[df_structure_summary['PDB'] == PDB_id]
        
        for model in structure:
            for chain in model:
                
                current_chain_df = current_PDB_df[current_PDB_df['pdbx_strand_id'].str.contains(chain.id)]
                row_count = len(current_chain_df)
                chain_id = chain.id
                
                if row_count == 1:
                    antibody_label = current_chain_df['antibody_label'].tolist()[0]
                    peptide_label = current_chain_df['peptide_label'].tolist()[0]

                    output_file_stem = "{}_{}_{}_{}".format(antibody_label, peptide_label, PDB_id, chain_id)
                    # print(FASTA_file_stem)
                    
                    if peptide_label == "heavy":
                        print("{} {} {} {}".format(PDB_id, chain_id, peptide_label, antibody_label))
                        sequence = ''
                        chain_type = 'H'
                        for polypeptide in ppb.build_peptides(chain):
                            sequence += polypeptide.get_sequence()
                            
                        FASTA_file_name = "{}.fasta".format(output_file_stem)
                        FASTA_file_path = Path(fasta_path, FASTA_file_name)

                        FASTA_file_path.write_text(">{}\n{}".format(chain_type, str(sequence)))
                        
                        # Make a dataframe to store all the renumber information
                        df_residue = pd.DataFrame()
                        
                        sequence_iter = ''
                        
                        for residue in chain:
                            threeletter = residue.get_resname()
                            
                            if threeletter not in residues_ignore_list:
                                hetero_field, resseq, icode = residue.get_id()
                                
                                try:
                                    oneletter = oneletter_dict[threeletter]
                                except:
                                    if threeletter not in exceptions_list:
                                        exceptions_list.append(threeletter)
                                    print(PDB_id, chain.get_id(), threeletter)

                                sequence_iter += oneletter

                                residue_dict = {
                                    'PDB': PDB_id,
                                    'chain': chain.get_id(),
                                    'threeletter': threeletter,
                                    'oneletter': oneletter,
                                    'resseq': resseq,
                                    'icode': icode
                                }

                                df_residue_dict = pd.DataFrame([residue_dict])

                                df_residue = pd.concat([df_residue, df_residue_dict], ignore_index=True)
                                
                        schemes = ['imgt', 'chothia', 'kabat']
                        for scheme in schemes:
                            
                            try:
                                abnumber_chain = Chain(sequence, scheme=scheme)
                            except Exception as e:
                                print(e)

                            abnumber_sequence = ''
                            abnumber_pos_list = []
                            abnumber_aa_list = []
                            new_numbers_list = []

                            for pos, aa in abnumber_chain:
                                # print(pos, aa)
                                abnumber_sequence += aa
                                abnumber_pos_list.append(pos)
                                abnumber_aa_list.append(aa)

                            # Sequences are the same

                            counter = 0
                            for index, row in df_residue.iterrows():
                                resseq_icode = str(row['resseq']) + str(row['icode'])
                                if counter < len(abnumber_aa_list):
                                    new_numbers_list.append(abnumber_pos_list[counter])
                                else:
                                    new_numbers_list.append('')
                                counter = counter + 1
                                # print(resseq_icode)

                            df_residue[scheme] = new_numbers_list

                
                
                        CSV_file_name = "{}.csv".format(output_file_stem)
                        CSV_output_file_path = Path(csv_path, 'renumbering_schemes', CSV_file_name)
                        df_residue.to_csv(CSV_output_file_path)                    
                    
                    if peptide_label == "light":
                        print("{} {} {} {}".format(PDB_id, chain_id, peptide_label, antibody_label))
                        sequence = ''
                        chain_type = 'L'
                        for polypeptide in ppb.build_peptides(chain):
                            sequence += polypeptide.get_sequence()

                        FASTA_file_name = "{}.fasta".format(output_file_stem)
                        FASTA_file_path = Path(fasta_path, FASTA_file_name)

                        FASTA_file_path.write_text(">{}\n{}".format(chain_type, str(sequence)))
                        
                        # Make a dataframe to store all the renumber information
                        df_residue = pd.DataFrame()
                        
                        sequence_iter = ''
                        
                        for residue in chain:
                            threeletter = residue.get_resname()
                            
                            if threeletter not in residues_ignore_list:
                                hetero_field, resseq, icode = residue.get_id()
                                
                                try:
                                    oneletter = oneletter_dict[threeletter]
                                except:
                                    if threeletter not in exceptions_list:
                                        exceptions_list.append(threeletter)
                                    print(PDB_id, chain.get_id(), threeletter)

                                sequence_iter += oneletter

                                residue_dict = {
                                    'PDB': PDB_id,
                                    'chain': chain.get_id(),
                                    'threeletter': threeletter,
                                    'oneletter': oneletter,
                                    'resseq': resseq,
                                    'icode': icode
                                }

                                df_residue_dict = pd.DataFrame([residue_dict])

                                df_residue = pd.concat([df_residue, df_residue_dict], ignore_index=True)
                                
                        schemes = ['imgt', 'chothia', 'kabat']
                        for scheme in schemes:
                            
                            abnumber_chain = ''

                            try:
                                abnumber_chain = Chain(sequence, scheme=scheme)
                            except Exception as e:
                                print(e)
                                print(residue_dict)

                            abnumber_sequence = ''
                            abnumber_pos_list = []
                            abnumber_aa_list = []
                            new_numbers_list = []

                            for pos, aa in abnumber_chain:
                                # print(pos, aa)
                                abnumber_sequence += aa
                                abnumber_pos_list.append(pos)
                                abnumber_aa_list.append(aa)

                            # Sequences are the same

                            counter = 0
                            for index, row in df_residue.iterrows():
                                resseq_icode = str(row['resseq']) + str(row['icode'])
                                if counter < len(abnumber_aa_list):
                                    new_numbers_list.append(abnumber_pos_list[counter])
                                else:
                                    new_numbers_list.append('')
                                counter = counter + 1
                                # print(resseq_icode)

                            df_residue[scheme] = new_numbers_list
                
                        CSV_file_name = "{}.csv".format(output_file_stem)
                        CSV_output_file_path = Path(csv_path, 'renumbering_schemes', CSV_file_name)
                        CSV_output_file_path.parents[0].mkdir(parents=True, exist_ok=True)
                        df_residue.to_csv(CSV_output_file_path)                
                    
                elif row_count > 1:
                    print('too many rows for {}{}'.format(PDB_id, chain_id))
                else:
                    # print('WARNING: 0 rows for {}{}'.format(PDB_id, chain_id))
                    continue

if __name__ == "__main__":
    main()