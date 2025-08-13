import pandas as pd
import mygene
import logging
import gprofiler

# Set up logging for better visibility
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_and_preprocess_data(uploaded_file):
    """
    Loads and pre-processes the CRISPR screen data from a TSV file.
    - Renames the first column to 'id' for consistency.
    - Explicitly converts key columns to numeric to prevent comparison errors.
    """
    try:
        data_df = pd.read_csv(uploaded_file, sep='\t')
        
        # Rename the first column to 'id'
        data_df.rename(columns={data_df.columns[0]: 'id'}, inplace=True)
        
        # Strip any leading/trailing whitespace from column names
        data_df.columns = data_df.columns.str.strip()
        
        # Explicitly convert LFC and FDR columns to numeric types
        numeric_cols = [
            'neg|lfc', 'pos|lfc', 
            'neg|fdr', 'pos|fdr', 
            'neg|score', 'pos|score',
            'num', 'neg|rank', 'pos|rank',
            'neg|goodsgrna', 'pos|goodsgrna'
        ]
        
        for col in numeric_cols:
            if col in data_df.columns:
                # Use errors='coerce' to turn non-numeric values into NaN
                data_df[col] = pd.to_numeric(data_df[col], errors='coerce')

        logging.info("CRISPR screen data loaded and preprocessed successfully.")
        return data_df
        
    except Exception as e:
        logging.error(f"Error loading and preprocessing CRISPR data: {e}", exc_info=True)
        return None

def merge_transcriptome_data(screen_df, transcriptome_file, screen_gene_col, transcriptome_gene_col, merge_type, transcriptome_has_header, organism):
    """
    Merges CRISPR screen data with transcriptome data.
    - Resolves gene aliases to ensure accurate merging.
    """
    try:
        # Load transcriptome data
        header_row = 0 if transcriptome_has_header else None
        transcriptome_df = pd.read_csv(transcriptome_file, sep=None, engine='python', header=header_row)
        
        # If no header, set column name for gene IDs
        if not transcriptome_has_header:
            transcriptome_df.rename(columns={transcriptome_df.columns[0]: transcriptome_gene_col}, inplace=True)
        
        # Clean up column names and gene IDs
        transcriptome_df.columns = [col.strip() for col in transcriptome_df.columns]
        transcriptome_df[transcriptome_gene_col] = transcriptome_df[transcriptome_gene_col].astype(str).str.strip()
        screen_df[screen_gene_col] = screen_df[screen_gene_col].astype(str).str.strip()

        # Resolve aliases for both datasets
        mg = mygene.MyGeneInfo()
        species_mapping = {'Mouse': 'mouse', 'Human': 'human', 'Rat': 'rat'}
        species_name = species_mapping.get(organism, 'mouse')
        
        # Resolve aliases for screen data
        screen_genes = screen_df[screen_gene_col].unique().tolist()
        screen_gene_info = mg.querymany(screen_genes, scopes='symbol,alias', species=species_name, as_dataframe=True)
        
        screen_gene_map = {}
        if not screen_gene_info.empty and 'query' in screen_gene_info.columns and 'symbol' in screen_gene_info.columns:
            screen_gene_map = {row['query']: row.get('symbol', row['query']) for _, row in screen_gene_info.iterrows()}
        screen_df['official_symbol'] = screen_df[screen_gene_col].map(screen_gene_map)
        screen_df['official_symbol'].fillna(screen_df[screen_gene_col], inplace=True)

        # Resolve aliases for transcriptome data
        transcriptome_genes = transcriptome_df[transcriptome_gene_col].unique().tolist()
        transcriptome_gene_info = mg.querymany(transcriptome_genes, scopes='symbol,alias', species=species_name, as_dataframe=True)

        transcriptome_gene_map = {}
        if not transcriptome_gene_info.empty and 'query' in transcriptome_gene_info.columns and 'symbol' in transcriptome_gene_info.columns:
            transcriptome_gene_map = {row['query']: row.get('symbol', row['query']) for _, row in transcriptome_gene_info.iterrows()}
        transcriptome_df['official_symbol'] = transcriptome_df[transcriptome_gene_col].map(transcriptome_gene_map)
        transcriptome_df['official_symbol'].fillna(transcriptome_df[transcriptome_gene_col], inplace=True)

        # Merge based on official symbols
        merged_df = pd.merge(
            screen_df, 
            transcriptome_df, 
            left_on='official_symbol', 
            right_on='official_symbol', 
            how=merge_type,
            suffixes=('_crispr', '_transcriptome')
        )
        
        logging.info("Data merged successfully.")
        return merged_df
    except Exception as e:
        logging.error(f"Error merging data: {e}", exc_info=True)
        return None

def filter_hits(df, selection_type, fdr_threshold, min_guides):
    """
    Filters the DataFrame based on selection type, FDR threshold, and minimum guides.
    This function relies on the `load_and_preprocess_data` function to ensure correct dtypes.
    """
    try:
        filtered_df = df.copy()

        # Check for existence of required columns
        required_cols = ['neg|fdr', 'pos|fdr', 'neg|goodsgrna', 'pos|goodsgrna']
        if not all(col in filtered_df.columns for col in required_cols):
            logging.error("Required columns for filtering (e.g., 'neg|fdr') are missing.")
            return pd.DataFrame()

        # Apply filtering logic
        if selection_type == 'negative':
            filtered_df = filtered_df[
                (filtered_df['neg|fdr'] < fdr_threshold) & 
                (filtered_df['neg|goodsgrna'] >= min_guides)
            ]
        elif selection_type == 'positive':
            filtered_df = filtered_df[
                (filtered_df['pos|fdr'] < fdr_threshold) & 
                (filtered_df['pos|goodsgrna'] >= min_guides)
            ]
        else: # 'both'
            filtered_df = filtered_df[
                ((filtered_df['neg|fdr'] < fdr_threshold) & 
                 (filtered_df['neg|goodsgrna'] >= min_guides)) |
                ((filtered_df['pos|fdr'] < fdr_threshold) & 
                 (filtered_df['pos|goodsgrna'] >= min_guides))
            ]
        
        # Sort by the minimum FDR
        if not filtered_df.empty:
            filtered_df['min_fdr'] = filtered_df[['neg|fdr', 'pos|fdr']].min(axis=1)
            filtered_df = filtered_df.sort_values(by='min_fdr')
            
        logging.info(f"Filtered to {len(filtered_df)} hits with '{selection_type}' selection.")
        return filtered_df
        
    except Exception as e:
        logging.error(f"Error filtering data: {e}", exc_info=True)
        return pd.DataFrame()

def perform_enrichment_analysis(gene_list, organism, p_value_threshold=0.05):
    """
    Performs enrichment analysis using g:Profiler.
    """
    try:
        gpro = gprofiler.GProfiler(return_dataframe=True)
        species_mapping = {'Mouse': 'mmusculus', 'Human': 'hsapiens', 'Rat': 'rnorvegicus'}
        species_name = species_mapping.get(organism, 'mmusculus')
        
        results = gpro.profile(
            organism=species_name,
            query=gene_list,
            no_evidences=False,
            user_threshold=p_value_threshold
        )
        
        if results.empty:
            logging.info("No significant enrichment terms found.")
            return {}

        results = results[results['p_value'] < p_value_threshold]
        
        # Split results into a dictionary of DataFrames by source
        grouped_results = {source: df for source, df in results.groupby('source')}
        
        logging.info("Enrichment analysis complete.")
        return grouped_results
        
    except Exception as e:
        logging.error(f"Error performing enrichment analysis: {e}", exc_info=True)
        return {}
        
def analyze_rnaseq_polarization(rnaseq_df, rnaseq_gene_col, significant_genes_list):
    """
    Analyzes RNA-seq data for significant CRISPR hits by polarization state.
    Reshapes the data and filters for the significant genes.
    """
    try:
        # Load and clean RNA-seq data
        rnaseq_df.columns = [col.strip() for col in rnaseq_df.columns]
        
        # Filter RNA-seq data to only include significant genes from CRISPR screen
        mg = mygene.MyGeneInfo()
        gene_info = mg.querymany(significant_genes_list, scopes='symbol,alias', species='mouse', as_dataframe=True)
        gene_map = {row['query']: row.get('symbol', row['query']) for _, row in gene_info.iterrows()}
        official_symbols = list(gene_map.values())
        
        filtered_rnaseq_df = rnaseq_df[rnaseq_df[rnaseq_gene_col].isin(official_symbols)].copy()
        
        # Add a column for official symbol
        filtered_rnaseq_df['official_symbol'] = filtered_rnaseq_df[rnaseq_gene_col].apply(lambda x: x if x in official_symbols else None)
        filtered_rnaseq_df.dropna(subset=['official_symbol'], inplace=True)
        
        # Melt the DataFrame to long format for easier plotting by polarization state
        # Get all columns except the gene column and any symbols added during merging
        value_vars = [col for col in filtered_rnaseq_df.columns if col not in [rnaseq_gene_col, 'official_symbol'] and '_' in col]
        
        melted_df = filtered_rnaseq_df.melt(
            id_vars=['official_symbol'], 
            value_vars=value_vars,
            var_name='Polarization_State', 
            value_name='Expression_Value'
        )
        
        # Convert Expression_Value to numeric, coercing errors
        melted_df['Expression_Value'] = pd.to_numeric(melted_df['Expression_Value'], errors='coerce').fillna(0)
        
        logging.info("RNA-seq polarization data analyzed successfully.")
        return melted_df

    except Exception as e:
        logging.error(f"Error analyzing RNA-seq polarization data: {e}", exc_info=True)
        return pd.DataFrame()

def analyze_rnaseq_data_timeseries(rnaseq_df, rnaseq_gene_col, rnaseq_time_col, rnaseq_value_col, significant_genes_list):
    """
    Analyzes RNA-seq data for significant CRISPR hits (time-series).
    """
    try:
        # Filter RNA-seq data to only include significant genes from CRISPR screen
        mg = mygene.MyGeneInfo()
        gene_info = mg.querymany(significant_genes_list, scopes='symbol,alias', species='mouse', as_dataframe=True)
        gene_map = {row['query']: row.get('symbol', row['query']) for _, row in gene_info.iterrows()}
        official_symbols = list(gene_map.values())
        
        filtered_rnaseq_df = rnaseq_df[rnaseq_df[rnaseq_gene_col].isin(official_symbols)].copy()
        
        # Add a column for official symbol
        filtered_rnaseq_df['official_symbol'] = filtered_rnaseq_df[rnaseq_gene_col].apply(lambda x: x if x in official_symbols else None)
        filtered_rnaseq_df.dropna(subset=['official_symbol'], inplace=True)
        
        # Convert time and value columns to numeric
        filtered_rnaseq_df[rnaseq_time_col] = pd.to_numeric(filtered_rnaseq_df[rnaseq_time_col], errors='coerce')
        filtered_rnaseq_df[rnaseq_value_col] = pd.to_numeric(filtered_rnaseq_df[rnaseq_value_col], errors='coerce')
        filtered_rnaseq_df.dropna(subset=[rnaseq_time_col, rnaseq_value_col], inplace=True)

        logging.info("RNA-seq time-series data analyzed successfully.")
        return filtered_rnaseq_df

    except Exception as e:
        logging.error(f"Error analyzing RNA-seq time-series data: {e}", exc_info=True)
        return pd.DataFrame()