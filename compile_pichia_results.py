import os
import zipfile
import pandas as pd
import h5py
import subprocess
import io
import tempfile

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(BASE_DIR, "results")
AGG_DIR = os.path.join(BASE_DIR, "aggregated_data")
OUTPUT_DIR = os.path.join(BASE_DIR, "aggregated_DEE_data")
ACCESSIONS_FILE = os.path.join(BASE_DIR, "kphaffii_accessions.txt")

# Organism info
ORG = "kphaffii"
SCIENTIFIC_NAME = "Komagataella phaffii"

def run_command(cmd):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running command: {cmd}\n{result.stderr}")
    return result.stdout

def get_metadata(accessions):
    print(f"Fetching metadata for organism '{SCIENTIFIC_NAME}' from ENA...")
    fields = (
        "run_accession,study_accession,experiment_accession,sample_accession,"
        "experiment_title,tax_id,scientific_name,library_strategy,library_source,"
        "library_selection,instrument_model,base_count,sample_alias,study_alias,"
        "secondary_study_accession,sample_description,study_title,description,"
        "library_construction_protocol,experimental_protocol"
    )
    # We fetch all transcriptomic data for this organism and then filter locally by our accessions
    url = f"https://www.ebi.ac.uk/ena/portal/api/search?query=scientific_name%3D%22{SCIENTIFIC_NAME.replace(' ', '%20')}%22%20AND%20library_source%3D%22TRANSCRIPTOMIC%22&result=read_run&fields={fields}&format=tsv&limit=0"
    cmd = f'curl -s "{url}"'
    output = run_command(cmd)
    if output.strip() and not output.startswith("Invalid"):
        try:
            df = pd.read_csv(io.StringIO(output), sep='\t')
            print(f"Total available metadata rows: {len(df)}")
            if 'run_accession' in df.columns:
                # Filter to our specific accessions
                df = df[df['run_accession'].isin(accessions)]
                return df
            else:
                print("Error: 'run_accession' column not found in ENA response.")
                print(f"Output starts with: {output[:100]}")
        except Exception as e:
            print(f"Error parsing ENA response: {e}")
            
    return pd.DataFrame()

def calculate_qc_summary(qc_data, num_genes=5393):
    # Replicate logic from dee_pipeline_functions.R
    try:
        # qc_data is a dict or Series
        num_reads_qc_pass = float(qc_data.get('NumReadsQcPass', 0))
        num_reads_qc_pass_per_gene = num_reads_qc_pass / num_genes
        
        qc_pass_rate_str = qc_data.get('QcPassRate', '0%').replace('%', '')
        qc_pass_rate = float(qc_pass_rate_str)
        
        star_uniq_map_rate_str = qc_data.get('STAR_UniqMapRate', '0%').replace('%', '')
        star_uniq_map_rate = float(star_uniq_map_rate_str)
        
        star_assign_rate_str = qc_data.get('STAR_AssignRate', '0%').replace('%', '')
        star_assign_rate = float(star_assign_rate_str)
        
        star_assigned_reads = float(qc_data.get('STAR_AssignedReads', 0))
        star_assigned_reads_per_gene = star_assigned_reads / num_genes
        
        kallisto_map_rate_str = qc_data.get('Kallisto_MapRate', '0%').replace('%', '')
        kallisto_map_rate = float(kallisto_map_rate_str)
        
        kallisto_mapped_reads = float(qc_data.get('Kallisto_MappedReads', 0))
        kallisto_mapped_reads_per_gene = kallisto_mapped_reads / num_genes
        
        fails = []
        warns = []
        
        # 1. NumReadsQcPass
        if num_reads_qc_pass_per_gene < 50: fails.append("1")
        elif num_reads_qc_pass_per_gene < 500: warns.append("1")
        
        # 2. QcPassRate
        if qc_pass_rate < 60: fails.append("2")
        elif qc_pass_rate < 80: warns.append("2")
        
        # 3. STAR_UniqMapRate
        if star_uniq_map_rate < 50: fails.append("3")
        elif star_uniq_map_rate < 70: warns.append("3")
        
        # 4. STAR_AssignRate
        if star_assign_rate < 40: fails.append("4")
        elif star_assign_rate < 60: warns.append("4")
        
        # 5. STAR_AssignedReads
        if star_assigned_reads_per_gene < 50: fails.append("5")
        elif star_assigned_reads_per_gene < 500: warns.append("5")
        
        # 6. Kallisto_MapRate
        if kallisto_map_rate < 40: fails.append("6")
        elif kallisto_map_rate < 60: warns.append("6")
        
        # 7. Kallisto_MappedReads
        if kallisto_mapped_reads_per_gene < 50: fails.append("7")
        elif kallisto_mapped_reads_per_gene < 500: warns.append("7")
        
        if fails:
            return f"FAIL({','.join(fails)})"
        if warns:
            return f"WARN({','.join(warns)})"
        return "PASS"
    except Exception as e:
        return f"ERROR({str(e)})"

def create_h5(filename, matrix_df, storage_mode='int'):
    print(f"Creating {filename}...")
    with h5py.File(filename, 'w') as f:
        data = matrix_df.values.T # (samples, genes)
        if storage_mode == 'int':
            dtype = 'int32'
        elif storage_mode == 'double':
            dtype = 'float64'
        else:
            # For QC, we'll store as fixed-length strings
            # Determine max length
            max_len = 0
            for val in data.flatten():
                max_len = max(max_len, len(str(val)))
            dtype = h5py.string_dtype(encoding='utf-8', length=max_len)
            data = data.astype(dtype)
            
        dset = f.create_dataset("bigmatrix", data=data, compression="gzip", compression_opts=9)
        
        # Rownames (Accessions)
        row_names = matrix_df.columns.values.astype(h5py.string_dtype())
        f.create_dataset("rownames", data=row_names)
        
        # Colnames (Gene/Tx IDs)
        col_names = matrix_df.index.values.astype(h5py.string_dtype())
        f.create_dataset("colnames", data=col_names)

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    with open(ACCESSIONS_FILE, 'r') as f:
        accessions = [line.strip() for line in f if line.strip()]
        
    print(f"Processing {len(accessions)} accessions...")
    
    gene_data_list = {}
    tx_data_list = {}
    qc_records = {}
    qc_summaries = {}
    
    available_accs = []

    for i, acc in enumerate(accessions):
        zip_path = os.path.join(RESULTS_DIR, f"{acc}.{ORG}.zip")
        if not os.path.exists(zip_path):
            zip_path = os.path.join(RESULTS_DIR, f"{acc}.zip")
            
        if os.path.exists(zip_path):
            try:
                with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                    # 1. Process Gene Counts (.se.tsv)
                    se_name = next((f for f in zip_ref.namelist() if f.endswith(f"{acc}.se.tsv")), None)
                    if se_name:
                        with zip_ref.open(se_name) as f:
                            # Skip 4 lines of STAR metadata
                            content = f.read().decode('utf-8').splitlines()
                            df_se = pd.read_csv(io.StringIO("\n".join(content[4:])), sep='\t', header=None, names=['ID', 'Count'], index_col=0)
                            gene_data_list[acc] = df_se['Count']
                    
                    # 2. Process Tx Counts (.ke.tsv)
                    ke_name = next((f for f in zip_ref.namelist() if f.endswith(f"{acc}.ke.tsv")), None)
                    if ke_name:
                        with zip_ref.open(ke_name) as f:
                            df_ke = pd.read_csv(f, sep='\t', index_col=0)
                            # Column 3 (0-indexed) is est_counts in the ke.tsv format: ID, length, eff_length, est_counts, tpm
                            tx_data_list[acc] = df_ke.iloc[:, 2] 
                    
                    # 3. Process QC (.qc)
                    qc_name = next((f for f in zip_ref.namelist() if f.endswith(f"{acc}.qc")), None)
                    if qc_name:
                        with zip_ref.open(qc_name) as f:
                            content = f.read().decode('utf-8')
                            qc_dict = {}
                            for line in content.strip().split('\n'):
                                if ':' in line:
                                    k, v = line.split(':', 1)
                                    qc_dict[k.strip()] = v.strip()
                            qc_records[acc] = qc_dict
                
                available_accs.append(acc)
                if i % 100 == 0 and i > 0: print(f"Processed {i} accessions...") # Progress hint
            except Exception as e:
                print(f"Error processing {acc} zip: {e}")
        else:
            # Try unzipped directory
            run_dir = os.path.join(RESULTS_DIR, acc)
            if os.path.exists(run_dir):
                try:
                    # Gene Counts
                    se_path = os.path.join(run_dir, f"{acc}.se.tsv")
                    if os.path.exists(se_path):
                        with open(se_path, 'r') as f:
                            content = f.readlines()
                            df_se = pd.read_csv(io.StringIO("".join(content[4:])), sep='\t', header=None, names=['ID', 'Count'], index_col=0)
                            gene_data_list[acc] = df_se['Count']
                    
                    # Tx Counts
                    ke_path = os.path.join(run_dir, f"{acc}.ke.tsv")
                    if os.path.exists(ke_path):
                        df_ke = pd.read_csv(ke_path, sep='\t', index_col=0)
                        tx_data_list[acc] = df_ke.iloc[:, 2]
                        
                    # QC
                    qc_path = os.path.join(run_dir, f"{acc}.qc")
                    if os.path.exists(qc_path):
                        with open(qc_path, 'r') as f:
                            qc_dict = {}
                            for line in f:
                                if ':' in line:
                                    k, v = line.split(':', 1)
                                    qc_dict[k.strip()] = v.strip()
                            qc_records[acc] = qc_dict
                    
                    available_accs.append(acc)
                except Exception as e:
                    print(f"Error processing {acc} directory: {e}")

    if not available_accs:
        print("No processed data found. Exiting.")
        return

    print(f"Creating matrices for {len(available_accs)} accessions...")
    gene_counts = pd.DataFrame(gene_data_list)
    tx_counts = pd.DataFrame(tx_data_list)
    
    # Calculate QC summaries
    for acc in available_accs:
        if acc in qc_records:
            qc_summaries[acc] = calculate_qc_summary(qc_records[acc], num_genes=len(gene_counts))

    # Create QC Matrix
    if qc_records:
        # Get all unique QC keys across all records
        first_acc = available_accs[0]
        all_keys = list(qc_records[first_acc].keys())
        qc_df = pd.DataFrame.from_dict(qc_records, orient='index')[all_keys].T
        create_h5(os.path.join(OUTPUT_DIR, "pichia_qc.h5"), qc_df, storage_mode='str')
    
    # Create SE and KE H5 files
    create_h5(os.path.join(OUTPUT_DIR, "pichia_se.h5"), gene_counts, storage_mode='int')
    create_h5(os.path.join(OUTPUT_DIR, "pichia_ke.h5"), tx_counts, storage_mode='double')
    
    # 4. Fetch metadata
    meta_df = get_metadata(available_accs)
    if not meta_df.empty:
        print(f"Metadata fetched for {len(meta_df)} accessions.")
        # Add QC_summary
        meta_df['QC_summary'] = meta_df['SRR_accession'].map(qc_summaries) if 'SRR_accession' in meta_df.columns else meta_df['run_accession'].map(qc_summaries)
        # Reorder columns to match E. coli format
        meta_df = meta_df.rename(columns={
            'run_accession': 'SRR_accession',
            'secondary_study_accession': 'SRP_accession', # Often SRP
            'experiment_accession': 'SRX_accession',
            'sample_accession': 'SRS_accession',
            'experiment_title': 'Experiment_title',
            'study_alias': 'GEO_series'
        })

        # If SRP_accession is still missing or NA, use study_accession
        if 'SRP_accession' in meta_df.columns:
            meta_df['SRP_accession'] = meta_df['SRP_accession'].fillna(meta_df['study_accession'])
        else:
            meta_df['SRP_accession'] = meta_df['study_accession']

        cols_to_front = [
            'SRR_accession', 'QC_summary', 'SRX_accession', 'SRS_accession', 'SRP_accession', 
            'Experiment_title', 'study_abstract', 'study_description', 'design_description', 
            'sample_description', 'sample_attribute', 'GEO_series'
        ]

        # Ensure all cols_to_front exist
        for col in cols_to_front:
            if col not in meta_df.columns:
                meta_df[col] = ""
                
        remaining_cols = [c for c in meta_df.columns if c not in cols_to_front]
        meta_df = meta_df[cols_to_front + remaining_cols]
        
        # Save metadata files
        print(f"Saving metadata to {OUTPUT_DIR}...")
        meta_df.to_csv(os.path.join(OUTPUT_DIR, "pichia_metadata.tsv"), sep='\t', index=False)
        meta_df[cols_to_front].to_csv(os.path.join(OUTPUT_DIR, "pichia_accessions.tsv"), sep='\t', index=False)
        meta_df[cols_to_front].to_csv(os.path.join(OUTPUT_DIR, "pichia_metadata.tsv.cut"), sep='\t', index=False)
        
        # Create pichia_srp.tsv
        print("Creating SRP metadata...")
        srp_accs = [s for s in meta_df['SRP_accession'].unique() if s and str(s) != 'nan']
        srp_list = []
        for srp in srp_accs:
            print(f"Fetching study metadata for {srp}...")
            # Use search API for study as well, it's more reliable
            url = f"https://www.ebi.ac.uk/ena/portal/api/search?query=study_accession%3D%22{srp}%22%20OR%20secondary_study_accession%3D%22{srp}%22&result=study&fields=study_accession,study_title,study_description,center_name,first_public&format=tsv&limit=1"
            cmd = f'curl -s "{url}"'
            output = run_command(cmd)
            if output.strip() and not output.startswith("Invalid"):
                try:
                    s_df = pd.read_csv(io.StringIO(output), sep='\t')
                    if not s_df.empty:
                        row = s_df.iloc[0]
                        srp_list.append({
                            'query_accession': srp,
                            'study_accession': row.get('study_accession', srp),
                            'title': row.get('study_title', ""),
                            'abstract': "",
                            'description': row.get('study_description', ""),
                            'study_type': "", # Not found in fields
                            'center': row.get('center_name', ""),
                            'submission_date': row.get('first_public', ""),
                            'error': "",
                            'GSE': "",
                            'URL': f"https://dee2.io/huge/pichia/{srp}_NA.zip"
                        })
                except Exception as e:
                    print(f"Error processing SRP {srp}: {e}")
        
        print(f"Total SRPs collected: {len(srp_list)}")
        if srp_list:
            srp_df = pd.DataFrame(srp_list)
            # Match the quoting style of ecoli_srp.tsv
            import csv
            srp_df.to_csv(os.path.join(OUTPUT_DIR, "pichia_srp.tsv"), sep='\t', index=False, quoting=csv.QUOTE_ALL)
    else:
        print("No metadata fetched. TSV files will not be created.")

    print(f"Compilation complete! Files saved in {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
