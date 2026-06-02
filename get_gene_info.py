import pandas as pd
import re

# Links collected from alife/bio/rnaseq/dee2/pipeline/volunteer_pipeline.sh
GENE_URLS = {
    "athaliana": "ftp://ftp.ensemblgenomes.org/pub/release-36/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.36.gtf.gz",
    "celegans": "ftp://ftp.ensembl.org/pub/release-90/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.90.gtf.gz",
    "dmelanogaster": "ftp://ftp.ensembl.org/pub/release-90/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.90.gtf.gz",
    "drerio": "ftp://ftp.ensembl.org/pub/release-90/gtf/danio_rerio/Danio_rerio.GRCz10.90.gtf.gz",
    "ecoli": "ftp://ftp.ensemblgenomes.org/pub/bacteria/release-36/gtf/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.36.gtf.gz",
    "hsapiens": "ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz",
    "mmusculus": "ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.gtf.gz",
    "rnorvegicus": "ftp://ftp.ensembl.org/pub/release-90/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.90.gtf.gz",
    "scerevisiae": "ftp://ftp.ensemblgenomes.org/pub/release-36/fungi/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.36.gtf.gz",
    "osativa": "ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/oryza_sativa/Oryza_sativa.IRGSP-1.0.59.gtf.gz",
    "zmays": "https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-59/plants/gtf/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.gtf.gz",
    "taestivum": "ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/triticum_aestivum/Triticum_aestivum.IWGSC.59.gtf.gz",
    "slycopersicum": "ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/solanum_lycopersicum/Solanum_lycopersicum.SL3.0.59.gtf.gz",
    "sbicolor": "ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/sorghum_bicolor/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.59.gtf.gz",
    "gmax": "ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/glycine_max/Glycine_max.Glycine_max_v2.1.59.gtf.gz",
    "ptrichocarpa": "ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/populus_trichocarpa/Populus_trichocarpa.Pop_tri_v4.59.gtf.gz",
    "vvinifera": "ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/vitis_vinifera/Vitis_vinifera.PN40024.v4.59.gtf.gz",
    "hvulgare": "ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/hordeum_vulgare/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.59.gtf.gz",
    "stuberosum": "ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/solanum_tuberosum/Solanum_tuberosum.SolTub_3.0.59.gtf.gz",
    "bdistachyon": "ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/brachypodium_distachyon/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.59.gtf.gz",
    "kphaffii": "https://ftp.ebi.ac.uk/ensemblgenomes/pub/fungi/release-62/gtf/komagataella_pastoris/Komagataella_pastoris.GCA_000027005.1.62.gtf.gz"
}

def extract_attribute(attr_str, key):
    """Extracts a value from the GTF attribute string."""
    match = re.search(f'{key} "([^"]+)"', attr_str)
    return match.group(1) if match else ""

def get_organism_data(organism: str):
    if organism not in GENE_URLS:
        print(f"Unknown organism: {organism}")
        return

    url = GENE_URLS[organism]
    print(f"Fetching {organism} data from {url}...")
    
    try:
        # Load GTF
        df_full = pd.read_csv(
            url, 
            sep="\t", 
            comment="#", 
            header=None,
            names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        )
        
        # 1. PROCESS GENE INFO
        df_genes = df_full[df_full['feature'] == 'gene'].copy()
        df_genes['Length'] = (df_genes['end'] - df_genes['start']) + 1
        df_genes['Gene_ID'] = df_genes['attribute'].apply(lambda x: extract_attribute(x, 'gene_id'))
        df_genes['Symbol'] = df_genes['attribute'].apply(lambda x: extract_attribute(x, 'gene_name') or extract_attribute(x, 'gene_id'))
        df_genes['Description'] = df_genes['attribute'].apply(lambda x: extract_attribute(x, 'description'))
        
        gene_info = df_genes[['Gene_ID', 'Symbol', 'Length', 'Description']].drop_duplicates()
        gene_file = f"{organism}_gene_info.csv"
        gene_info.to_csv(gene_file, index=False)
        print(f"Saved {gene_file} ({len(gene_info)} genes)")

        # 2. PROCESS TRANSCRIPT INFO
        # Look for 'transcript' or 'mRNA' features
        df_trans = df_full[df_full['feature'].isin(['transcript', 'mRNA'])].copy()
        
        if len(df_trans) == 0:
            # Some bacterial GTFs might only have CDS/gene.
            print(f"Warning: No transcript features found for {organism}. Checking CDS as proxy...")
            df_trans = df_full[df_full['feature'] == 'CDS'].copy()

        df_trans['Transcript_ID'] = df_trans['attribute'].apply(lambda x: extract_attribute(x, 'transcript_id') or extract_attribute(x, 'protein_id'))
        df_trans['Gene_ID'] = df_trans['attribute'].apply(lambda x: extract_attribute(x, 'gene_id'))
        df_trans['Symbol'] = df_trans['attribute'].apply(lambda x: extract_attribute(x, 'gene_name'))
        
        # Fallback to Gene_ID if Symbol is empty
        df_trans.loc[df_trans['Symbol'] == "", 'Symbol'] = df_trans['Gene_ID']
        
        # For transcripts, 'Length' in GTF is genomic span.
        df_trans['Length'] = (df_trans['end'] - df_trans['start']) + 1
        
        trans_info = df_trans[['Transcript_ID', 'Gene_ID', 'Symbol', 'Length']].drop_duplicates()
        # Remove entries with empty IDs
        trans_info = trans_info[trans_info['Transcript_ID'] != ""]
        
        trans_file = f"{organism}_transcript_info.csv"
        trans_info.to_csv(trans_file, index=False)
        print(f"Saved {trans_file} ({len(trans_info)} transcripts)")
        
    except Exception as e:
        print(f"Error processing {organism}: {e}")

# Fetch primary organisms
for org in ["kphaffii", "ecoli", "scerevisiae"]:
    get_organism_data(org)
