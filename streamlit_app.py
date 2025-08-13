import streamlit as st
import pandas as pd
import data_analysis # The updated data_analysis.py module
import google.generativeai as genai
import logging
import time
import plotly.express as px

# Set up logging for better visibility
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- Streamlit App Configuration ---
st.set_page_config(
    page_title="CRISPR Screen Analysis & AI Report",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("CRISPR Screen Data Analysis & AI Report Generator")


# --- CONSOLIDATED SESSION STATE INITIALIZATION ---
# This block ensures all variables exist before the app tries to use them.
if 'processed_data' not in st.session_state:
    st.session_state['processed_data'] = None
if 'filtered_hits' not in st.session_state:
    st.session_state['filtered_hits'] = pd.DataFrame()
if 'main_enrichment_results' not in st.session_state:
    st.session_state['main_enrichment_results'] = {}
if 'ai_report_main' not in st.session_state:
    st.session_state['ai_report_main'] = ""
if 'ai_report_targeted' not in st.session_state:
    st.session_state['ai_report_targeted'] = {}
if 'gene_list_for_enrichment' not in st.session_state:
    st.session_state['gene_list_for_enrichment'] = []
if 'targeted_pathway_genes' not in st.session_state:
    st.session_state['targeted_pathway_genes'] = []
if 'rnaseq_analyzed_hits_timeseries' not in st.session_state:
    st.session_state['rnaseq_analyzed_hits_timeseries'] = pd.DataFrame()
if 'rnaseq_analyzed_hits_polarization' not in st.session_state:
    st.session_state['rnaseq_analyzed_hits_polarization'] = pd.DataFrame()
if 'rnaseq_report' not in st.session_state:
    st.session_state['rnaseq_report'] = ""

# AI report context variables with default values
if 'user_research_focus' not in st.session_state:
    st.session_state['user_research_focus'] = "How do N-linked glycosylation pathways and specific N-glycosylated proteins influence the molecular mechanisms of ADCP in macrophages?"
if 'user_experimental_system' not in st.session_state:
    st.session_state['user_experimental_system'] = "genome-scale CRISPR-Cas9 knockout screen in primary murine bone marrow-derived macrophages (BMDMs) identifying genetic regulators of Antibody-Dependent Cellular Phagocytosis (ADCP) when macrophages encounter opsonized red blood cells (RBCs)."
if 'user_phenotype_interpretation' not in st.session_state:
    st.session_state['user_phenotype_interpretation'] = "A negative Log2 Fold Change (LFC) indicates that the knockout (loss of function) of a particular gene leads to an enhancement of the measured phenotype (ADCP), whereas a positive LFC suggests that gene knockout suppresses the phenotype."
if 'user_specific_areas' not in st.session_state:
    st.session_state['user_specific_areas'] = "N-linked Glycan Biosynthesis, Processing, and Quality Control, Trafficking and Presentation of N-glycosylated Proteins, Key N-glycosylated Proteins in Immune Recognition and ADCP Co-regulation, Lectins and Glycan-binding Proteins (GBPs) in Macrophages"
if 'user_priority_pathways' not in st.session_state:
    st.session_state['user_priority_pathways'] = ""


# --- Sidebar: Data Upload ---
st.sidebar.header("Upload Data")
uploaded_crispr_file = st.sidebar.file_uploader(
    "Upload CRISPR Screen Data (TSV, no header)", type=["tsv"]
)
organism_selection = st.sidebar.selectbox(
    "Select Organism:",
    ('Mouse', 'Human', 'Rat'),
    index=0 # Default to Mouse
)

if uploaded_crispr_file:
    with st.spinner("Loading and preprocessing CRISPR data..."):
        st.session_state['processed_data'] = data_analysis.load_and_preprocess_data(uploaded_crispr_file)
        if st.session_state['processed_data'] is not None:
            st.sidebar.success("CRISPR Screen Data loaded!")
        else:
            st.sidebar.error("Failed to load CRISPR data. Please check file format.")
            st.session_state['processed_data'] = None

# --- Sidebar: Merge Transcriptome Data (Optional) ---
st.sidebar.header("Merge Transcriptome Data (Optional)")
uploaded_transcriptome_file = st.sidebar.file_uploader(
    "Upload Transcriptome Data (CSV/TSV/TXT)", type=["csv", "tsv", "txt"]
)
transcriptome_has_header = st.sidebar.checkbox("Transcriptome file has header?", value=True)
transcriptome_gene_id_col = st.sidebar.text_input(
    "Transcriptome Gene ID Column Name:", value="Genenames",
    help="The column in your transcriptome file containing gene identifiers for merging (e.g., 'Genenames', 'gene_symbol')."
)
merge_option = st.sidebar.selectbox(
    "Merge Type:",
    options=['left', 'inner', 'outer', 'right'],
    index=0, # 'left' as default
    help="How to combine data:\n"
         "- Left (default): Keep all CRISPR genes, add matching transcriptome data.\n"
         "- Inner: Keep only genes present in BOTH datasets.\n"
         "- Outer: Keep all genes from either dataset.\n"
         "- Right: Keep all transcriptome genes, add matching CRISPR data."
)

if uploaded_transcriptome_file and st.session_state['processed_data'] is not None:
    if st.sidebar.button("Perform Data Merge", key="merge_button"):
        with st.spinner("Merging transcriptome data..."):
            merged_df = data_analysis.merge_transcriptome_data(
                screen_df=st.session_state['processed_data'],
                transcriptome_file=uploaded_transcriptome_file,
                screen_gene_col='id',
                transcriptome_gene_col=transcriptome_gene_id_col,
                merge_type=merge_option,
                transcriptome_has_header=transcriptome_has_header,
                organism=organism_selection
            )
            if merged_df is not None:
                st.session_state['processed_data'] = merged_df
                st.sidebar.success("Transcriptome data merged successfully!")
            else:
                st.sidebar.error("Failed to merge transcriptome data.")
elif uploaded_transcriptome_file and st.session_state['processed_data'] is None:
    st.sidebar.warning("Upload CRISPR data first to enable transcriptome merging.")


# --- Sidebar: Filters ---
st.sidebar.header("Filters")
st.sidebar.write("Adjust filters and click 'Run Analysis' to process data.")
hit_type = st.sidebar.selectbox(
    "Hit Type:",
    ('negative', 'positive', 'both'),
    help="Select to analyze negative selection hits, positive selection hits, or both."
)
fdr_threshold = st.sidebar.number_input(
    "FDR Threshold:",
    min_value=0.0,
    max_value=1.0,
    value=0.05,
    step=0.01,
    format="%.3f",
    help="False Discovery Rate (FDR) threshold for statistical significance."
)
min_guides = st.sidebar.number_input(
    "Minimum Guides Per Gene:",
    min_value=1,
    max_value=10,
    value=4,
    step=1,
    help="Minimum number of guides required for a gene to be considered a hit."
)
p_value_enrichment = st.sidebar.number_input(
    "Enrichment P-value Threshold:",
    min_value=0.0,
    max_value=1.0,
    value=0.05,
    step=0.01,
    format="%.3f",
    key="enrichment_p_value_input",
    help="P-value threshold for the gProfiler enrichment analysis."
)
st.sidebar.markdown("---")


# --- Main Page: Data & Analysis Results ---

if st.session_state['processed_data'] is not None:
    st.header("Raw Preprocessed Data")
    st.write("First 5 rows of the loaded CRISPR screen data:")
    st.dataframe(st.session_state['processed_data'].head())
    st.markdown("---")

    if st.button("Run Analysis", key="run_analysis_button"):
        with st.spinner("Filtering hits and performing gene enrichment analysis..."):
            if 'id' not in st.session_state['processed_data'].columns:
                st.error("Error: The loaded CRISPR data does not have a column named 'id'. "
                         "Please check your input TSV file.")
            else:
                # Step 1: Filter hits
                filtered_df = data_analysis.filter_hits(
                    st.session_state['processed_data'],
                    selection_type=hit_type,
                    fdr_threshold=fdr_threshold,
                    min_guides=min_guides
                )
                st.session_state['filtered_hits'] = filtered_df

                st.subheader("Filtered CRISPR Hits")
                if not st.session_state['filtered_hits'].empty:
                    st.session_state['gene_list_for_enrichment'] = st.session_state['filtered_hits']['id'].tolist()
                    st.success(f"Found {len(st.session_state['filtered_hits'])} significant hits.")
                    st.dataframe(st.session_state['filtered_hits'])
                    
                    # Step 2: Perform enrichment analysis
                    st.markdown("---")
                    st.subheader("Main Gene Enrichment Analysis")
                    st.info("The tables below show the top enriched pathways and terms.")
                    st.session_state['main_enrichment_results'] = data_analysis.perform_enrichment_analysis(
                        gene_list=st.session_state['gene_list_for_enrichment'],
                        organism=organism_selection,
                        p_value_threshold=p_value_enrichment
                    )
                    
                    if st.session_state['main_enrichment_results']:
                        tab_names = list(st.session_state['main_enrichment_results'].keys())
                        st.session_state['targeted_pathway_genes'] = [
                            term for source in st.session_state['main_enrichment_results'].values()
                            for term in source['name'].tolist()
                        ]
                        tabs = st.tabs(tab_names)
                        for i, tab in enumerate(tabs):
                            with tab:
                                source_name = tab_names[i]
                                df_results = st.session_state['main_enrichment_results'][source_name]
                                if not df_results.empty:
                                    st.subheader(f"Results for {source_name}")
                                    st.dataframe(df_results)
                                else:
                                    st.warning(f"No significant enrichment found for {source_name} with the current settings.")
                    else:
                        st.warning("No significant enrichment terms found. Try adjusting the p-value threshold or filters.")
                else:
                    st.session_state['gene_list_for_enrichment'] = []
                    st.session_state['main_enrichment_results'] = {}
                    st.warning("No significant hits found with the current filters. Please adjust your filters and try again.")
    
    # Re-display results if they already exist in session state on app reload
    if not st.session_state['filtered_hits'].empty:
        st.subheader("Filtered CRISPR Hits")
        st.success(f"Found {len(st.session_state['filtered_hits'])} significant hits.")
        st.dataframe(st.session_state['filtered_hits'])
        st.markdown("---")

    if st.session_state['main_enrichment_results']:
        st.subheader("Main Gene Enrichment Analysis")
        tab_names = list(st.session_state['main_enrichment_results'].keys())
        tabs = st.tabs(tab_names)
        for i, tab in enumerate(tabs):
            with tab:
                source_name = tab_names[i]
                df_results = st.session_state['main_enrichment_results'][source_name]
                st.subheader(f"Results for {source_name}")
                st.dataframe(df_results)
    
    st.markdown("---")


# --- Main Page: AI Report Prompts ---
st.header("Report Question Prompts")
st.write("This information will be used by the AI to generate a more specific and relevant research report.")
st.session_state['user_research_focus'] = st.text_input(
      "What is the core research question or focus of this study?",
      st.session_state['user_research_focus'],
      key="user_research_focus_input"
)
st.session_state['user_experimental_system'] = st.text_input(
      "Describe your experimental system (e.g., 'CRISPR screen in human T-cells'):",
      st.session_state['user_experimental_system'],
      key="user_experimental_system_input"
)
st.session_state['user_phenotype_interpretation'] = st.text_area(
      "How should Log2 Fold Change (LFC) be interpreted for your screen?",
      st.session_state['user_phenotype_interpretation'],
      key="user_phenotype_interpretation_input"
)
st.session_state['user_specific_areas'] = st.text_area(
    "Specific Areas of Interest for Gene Identification and Mechanistic Insight:",
    st.session_state['user_specific_areas'],
    key="user_specific_areas_input"
)
st.session_state['user_priority_pathways'] = st.text_area(
    "Specific Pathways for Detailed Analysis:",
    st.session_state['user_priority_pathways'],
    key="user_priority_pathways_input"
)
st.markdown("---")


# --- Main Page: Main AI Report ---
st.header("AI Report")
if not st.session_state['filtered_hits'].empty and st.session_state['main_enrichment_results']:
    if st.button("Generate Main AI Report"):
        if "GOOGLE_API_KEY" not in st.secrets:
            st.error("Google API key not found in .streamlit/secrets.toml. Please add it to generate the AI report.")
        else:
            with st.spinner("Generating AI report... This may take a few minutes."):
                try:
                    genai.configure(api_key=st.secrets["GOOGLE_API_KEY"])
                    model = genai.GenerativeModel('gemini-2.5-pro')

                    filtered_gene_data_for_prompt = st.session_state['filtered_hits'][['id', 'neg|lfc', 'pos|lfc', 'neg|fdr', 'pos|fdr']].to_markdown(index=False)
                    final_enrichment_summary_text = ""
                    for source, df in st.session_state['main_enrichment_results'].items():
                        final_enrichment_summary_text += f"### {source}\n"
                        if not df.empty:
                            final_enrichment_summary_text += df.head(10).to_markdown(index=False) + "\n\n"
                        else:
                            final_enrichment_summary_text += "No significant terms.\n\n"

                    report_fdr_threshold = fdr_threshold
                    
                    ai_prompt = f"""
                    You are an expert bioinformatician and molecular biologist, tasked with generating a comprehensive research report.
                    The report will analyze gene enrichment data to shed light on **{st.session_state['user_research_focus']}**. The data originates from a **{st.session_state['user_experimental_system']}**.

                    **Core Research Question for this Report:**
                    {st.session_state['user_research_focus']}

                    **Specific Areas of Interest for Gene Identification and Mechanistic Insight:**
                    {st.session_state['user_specific_areas']}

                    **Interpretation of Log2 Fold Change (LFC):**
                    {st.session_state['user_phenotype_interpretation']}

                    **Provided Gene Enrichment Data:**
                    (This data represents significant pathways/GO terms from your screen analysis. Please focus your interpretation primarily on these provided results. If you discuss other genes/pathways for contextual relevance, explicitly state they were not a significant hit in *this specific screen data*.)

                    ```
                    {final_enrichment_summary_text}
                    ```

                    **All Significant Genes from Screen with LFC and FDR:**
                    (This data lists all genes that passed your initial filtering. Use this to populate the LFC and FDR columns in Section IV, and to select relevant ones for discussion.)

                    ```
                    {filtered_gene_data_for_prompt}
                    ```

                    **Specific Pathways for Detailed Analysis:**
                    {st.session_state['user_priority_pathways'] if st.session_state['user_priority_pathways'] else 'No specific pathways were prioritized by the user. The analysis will focus on the most statistically significant terms.'}
                    If any of the "Specific Pathways for Detailed Analysis" are present in the provided enrichment data, ensure they are explicitly discussed and interpreted in the relevant sections below, even if they are not among the top-ranked terms. Relate them directly to the Core Research Question and mention their significance (e.g., FDR).

                    **Desired Report Structure and Analysis:**
                    The report should go beyond a simple list of enriched genes and provide a mechanistic interpretation.

                    I. Executive Summary
                    * Provide a brief overview of key findings and their implications for the study's core question.
                    * Highlight central findings related to the user's specific areas of interest, mentioning categories of genes and their LFC patterns.
                    * Discuss broader implications, e.g., how disruptions might influence cellular signals or processes.

                    II. Methodology
                    * Describe the origin of the screen data (based on user's experimental system description).
                    * Explain the meaning of Log2 Fold Change (LFC) and False Discovery Rate (FDR) in this context, using the user's interpretation.
                    * Mention the significance threshold used (FDR < {report_fdr_threshold}).
                    * Briefly describe the bioinformatics tools and databases used for analysis (e.g., gProfiler for enrichment, and conceptually KEGG, Reactome, Gene Ontology as underlying knowledge bases, and STRING for protein-protein interaction insights). Do not include numbered citations.

                    III. Overall Pathway Enrichment
                    * Discuss the overall influence of enriched pathways and related cellular processes based on the provided data in "Provided Gene Enrichment Data".
                    * If distinct conditions were provided in the "Provided Gene Enrichment Data" section, perform a comparative examination, noting quantitative differences and implications. Otherwise, focus on the data provided.
                    * Highlight consistent patterns in LFC values (referring to the "All Significant Genes" data).
                    * Interpret broader implications related to cellular health or fundamental processes.
                    * Present a summary table (Markdown format) of the top enriched pathways and GO terms from the most relevant dataset (from "Provided Gene Enrichment Data").
                        * Table Columns: "Pathway/GO Term Category", "Pathway/GO Term Name", "Representative Genes (FDR < {report_fdr_threshold})"
                        * For "Representative Genes", select 3-7 most relevant genes from *the provided enrichment data* for that term, prioritizing those aligned with the user's specific areas of interest.

                    IV. Targeted Gene Analysis: Specific Regulators
                    This section should provide a detailed examination of specific candidate genes identified from the screen (using the "All Significant Genes" data, FDR < {report_fdr_threshold}), categorizing them according to the user's "Specific Areas of Interest."
                    For each gene discussed:
                        * Include its LFC and FDR (from the "All Significant Genes" data, clearly indicating which condition if multiple are present).
                        * Describe its known function and proposed mechanistic link to the study's focus based on general biological knowledge.
                        * Highlight genes with strong direct or compelling indirect evidence.
                        * If you mention a gene that was *not* significant in the provided screen data for contextual relevance, explicitly state that it was not a significant hit in *this screen data*.

                    V. Protein-Protein Interaction Networks (Conceptual Discussion)
                    * Discuss conceptual interactions among highly enriched genes, drawing on the provided enrichment data and the list of significant genes.
                    * Describe what a network diagram (e.g., from STRING) would likely reveal about functional modules or connections, without generating an actual image.

                    VI. Hypothesis Generation and Future Directions
                    * Propose novel, testable hypotheses regarding the precise roles of identified pathways or specific genes in the study's focus, leveraging the provided data and your biological knowledge.
                    * Suggest potential experimental validation strategies.

                    VII. Discussion of Limitations
                    * Acknowledge any limitations of the analysis (e.g., reliance on existing databases, specificity of the initial enrichment data, potential for off-target effects in CRISPR screen, inability to provide real-time citations for specific papers, general limitations of LLM analysis).

                    **Tone:** Academic, analytical, insightful, and comprehensive.

                    Begin the report now:
                    """
                    response = model.generate_content(ai_prompt)
                    st.session_state['ai_report_main'] = response.text
                    st.success("Main AI report generated!")
                except Exception as e:
                    st.session_state['ai_report_main'] = f"Error generating AI report: {e}. Please try again or simplify the prompt/data."
                    st.error("Error generating AI report. Check logs for details.")
                    logging.error(f"Gemini API Error: {e}", exc_info=True)
            time.sleep(1)

    if st.session_state['ai_report_main']:
        st.subheader("Generated Main Research Report")
        st.markdown(st.session_state['ai_report_main'])
        st.download_button(
            label="Download Main AI Report",
            data=st.session_state['ai_report_main'].encode('utf-8'),
            file_name="CRISPR_main_AI_report.md",
            mime="text/markdown",
        )
else:
    st.info("Upload data, adjust filters, and click 'Run Analysis' to generate a main AI report.")

st.markdown("---")


# --- Main Page: Targeted Pathway Analysis (Optional) ---
st.header("Targeted Pathway Analysis (Optional)")
if st.session_state['main_enrichment_results'] and not st.session_state['filtered_hits'].empty:
    st.write("Generate a detailed AI report for specific pathways of interest.")
    pathway_selection_method = st.radio(
        "Select pathway by:",
        ('Identified Pathways', 'Manual Input'),
        key='pathway_selection_method'
    )
    selected_pathway_name = ""
    if pathway_selection_method == 'Identified Pathways':
        if st.session_state['targeted_pathway_genes']:
            selected_pathway_name = st.selectbox(
                "Choose a specific pathway name from enriched results:",
                options=[""] + sorted(list(set(st.session_state['targeted_pathway_genes']))),
                index=0,
                key="select_pathway_from_list"
            )
        else:
            st.warning("No enriched pathways identified yet. Run main analysis first.")
    else:
        selected_pathway_name = st.text_input(
            "Enter specific pathway name (e.g., 'immune response', 'MAPK signaling'):",
            key="manual_pathway_input"
        )
    
    st.sidebar.subheader("Targeted Pathway Filters")
    pathway_fdr_threshold = st.sidebar.number_input(
        "Pathway FDR Threshold (for specific terms):",
        min_value=0.0, max_value=1.0, value=0.05, step=0.01, format="%.3f",
        key="pathway_fdr_input",
        help="FDR threshold for the individual pathway term itself within the gProfiler results."
    )
    pathway_p_value_threshold = st.sidebar.number_input(
        "Pathway P-value Threshold (for specific terms):",
        min_value=0.0, max_value=1.0, value=0.05, step=0.01, format="%.3f",
        key="pathway_pvalue_input",
        help="P-value threshold for the individual pathway term itself within the gProfiler results."
    )
    min_genes_in_pathway = st.sidebar.number_input(
        "Min Genes in Pathway Intersection:",
        min_value=1, value=2, step=1,
        key="min_genes_pathway_input",
        help="Minimum number of your significant genes that must intersect with the pathway's genes."
    )

    if st.button("Generate AI Report for Targeted Pathway") and selected_pathway_name:
        if "GOOGLE_API_KEY" not in st.secrets:
            st.error("Google API key not found in .streamlit/secrets.toml. Please add it to generate the AI report.")
        else:
            with st.spinner(f"Generating AI report for {selected_pathway_name}..."):
                try:
                    genai.configure(api_key=st.secrets["GOOGLE_API_KEY"])
                    model = genai.GenerativeModel('gemini-2.5-pro')
                    pathway_details = None
                    genes_in_pathway_from_gprofiler = []
                    for source, df_results in st.session_state['main_enrichment_results'].items():
                        if not df_results.empty:
                            matched_pathways = df_results[
                                (df_results['name'].str.contains(selected_pathway_name, case=False, na=False)) &
                                (df_results['fdr'] < pathway_fdr_threshold) &
                                (df_results['p_value'] < pathway_p_value_threshold) &
                                (df_results['intersection_size'] >= min_genes_in_pathway)
                            ]
                            if not matched_pathways.empty:
                                pathway_details = matched_pathways.iloc[0].to_dict()
                                genes_in_pathway_from_gprofiler = [g.strip() for g in pathway_details.get('genes', '').split(',') if g.strip()]
                                logging.info(f"Found pathway details for '{selected_pathway_name}' in source '{source}'.")
                                break
                    if pathway_details:
                        targeted_report_content = ""
                        user_context_str = ""
                        if st.session_state['user_research_focus']:
                            user_context_str += f"**Research Focus:** {st.session_state['user_research_focus']}\n"
                        if st.session_state['user_experimental_system']:
                            user_context_str += f"**Experimental System:** {st.session_state['user_experimental_system']}\n"
                        if st.session_state['user_phenotype_interpretation']:
                            user_context_str += f"**Phenotype Interpretation:** {st.session_state['user_phenotype_interpretation']}\n\n"
                        targeted_report_content += user_context_str
                        targeted_report_content += f"Targeted Pathway: {selected_pathway_name}\n\n"
                        targeted_report_content += "Pathway Details from gProfiler:\n"
                        for k, v in pathway_details.items():
                            targeted_report_content += f"- {k}: {v}\n"
                        targeted_report_content += "\nOverall CRISPR Screen Significant Hits (Subset intersecting this pathway):\n"
                        significant_genes_in_this_pathway_df = st.session_state['filtered_hits'][
                            st.session_state['filtered_hits']['id'].isin(genes_in_pathway_from_gprofiler)
                        ]
                        if not significant_genes_in_this_pathway_df.empty:
                            targeted_report_content += significant_genes_in_this_pathway_df.to_markdown(index=False) + "\n\n"
                        else:
                            targeted_report_content += "No significant genes from your screen directly intersect with this pathway based on current filters and gProfiler intersection.\n\n"
                        if st.session_state['processed_data'] is not None and len(st.session_state['processed_data'].columns) > 14:
                            merged_genes_in_pathway_with_transcriptome_df = st.session_state['processed_data'][
                                st.session_state['processed_data']['id'].isin(significant_genes_in_this_pathway_df['id'].tolist())
                            ]
                            if not merged_genes_in_pathway_with_transcriptome_df.empty:
                                targeted_report_content += "Expression/Transcriptome Data for intersecting significant genes:\n"
                                transcriptome_cols_for_report = [col for col in merged_genes_in_pathway_with_transcriptome_df.columns if '_transcriptome' in col or col == transcriptome_gene_id_col]
                                common_transcriptome_cols = ['aveFLM', 'aveFLM-GMT']
                                for col in common_transcriptome_cols:
                                    if col in merged_genes_in_pathway_with_transcriptome_df.columns and col not in transcriptome_cols_for_report:
                                        transcriptome_cols_for_report.append(col)
                                if transcriptome_cols_for_report:
                                    targeted_report_content += merged_genes_in_pathway_with_transcriptome_df[transcriptome_cols_for_report].head(10).to_markdown(index=False) + "\n\n"
                                else:
                                    targeted_report_content += "No relevant transcriptome columns found for these genes.\n\n"
                            else:
                                targeted_report_content += "No transcriptome data available for the identified significant genes in this pathway.\n\n"

                        custom_targeted_prompt = st.text_area(
                            "Custom prompt for Targeted AI Report:",
                            value=f"Provide a detailed biological interpretation of the pathway '{selected_pathway_name}' based on the provided CRISPR screen data, its enrichment details, and specific intersecting significant genes. Explain its known biological roles, connection to the screen's context (e.g., macrophage polarization if relevant), and implications of the gene hits within this pathway. Include relevant expression data if available. Focus on key genes and their roles.",
                            height=150,
                            key="targeted_prompt_input"
                        )
                        response = model.generate_content(
                            f"{custom_targeted_prompt}\n\nData to analyze:\n{targeted_report_content}"
                        )
                        st.session_state['ai_report_targeted'][selected_pathway_name] = response.text
                        st.success(f"AI report generated for '{selected_pathway_name}'!")
                    else:
                        st.warning(f"Could not find details for pathway '{selected_pathway_name}' with current filters. Adjust pathway name, FDR, P-value, or min genes. Ensure it was found in the main enrichment step.")
                except Exception as e:
                    st.session_state['ai_report_targeted'][selected_pathway_name] = f"Error generating AI report for '{selected_pathway_name}': {e}. Please try again or simplify the prompt/data."
                    st.error(f"Error generating AI report for '{selected_pathway_name}'. Check logs for details.")
                    logging.error(f"Gemini API Error for targeted report: {e}", exc_info=True)
            time.sleep(1)

    if st.session_state['ai_report_targeted']:
        st.subheader("Generated Targeted Pathway Reports")
        for pathway, report_text in st.session_state['ai_report_targeted'].items():
            with st.expander(f"Report for: {pathway}"):
                st.markdown(report_text)
                st.download_button(
                    label=f"Download Report for {pathway}",
                    data=report_text.encode('utf-8'),
                    file_name=f"{pathway.replace(' ', '_')}_AI_report.md",
                    mime="text/markdown",
                    key=f"download_targeted_{pathway}"
                )
else:
    st.info("Run main enrichment analysis to enable targeted pathway analysis.")

st.markdown("---")


# --- Main Page: Optional RNA-seq Data Analysis ---
st.header("Optional RNA-seq Data Analysis")
st.write("Upload RNA-seq data to analyze expression trends for your significant CRISPR hits.")
uploaded_rnaseq_file = st.file_uploader(
    "Upload RNA-seq Data (CSV/TSV)",
    type=["csv", "tsv"],
    key="rnaseq_uploader"
)

if uploaded_rnaseq_file and not st.session_state['filtered_hits'].empty:
    st.session_state['rna_seq_data'] = pd.read_csv(uploaded_rnaseq_file, sep='\t' if uploaded_rnaseq_file.name.endswith('.tsv') else ',')
    rnaseq_cols = st.session_state['rna_seq_data'].columns.tolist()

    analysis_type = st.radio(
        "Select RNA-seq Analysis Type:",
        ('Time-series Analysis', 'Polarization State Comparison'),
        help="Choose to compare expression over time or between different conditions/polarization states."
    )

    if analysis_type == 'Time-series Analysis':
        rnaseq_gene_col = st.selectbox("Gene Name Column:", rnaseq_cols, index=0, key="rnaseq_gene_col")
        rnaseq_time_col = st.selectbox("Time Point Column:", rnaseq_cols, key="rnaseq_time_col")
        rnaseq_value_col = st.selectbox("Expression Value Column:", rnaseq_cols, key="rnaseq_value_col")

        if st.button("Run Time-series Analysis", key="run_timeseries_btn"):
            with st.spinner("Analyzing RNA-seq time-series data..."):
                st.session_state['rnaseq_analyzed_hits_timeseries'] = data_analysis.analyze_rnaseq_data_timeseries(
                    st.session_state['rna_seq_data'],
                    rnaseq_gene_col,
                    rnaseq_time_col,
                    rnaseq_value_col,
                    st.session_state['gene_list_for_enrichment']
                )

                if not st.session_state['rnaseq_analyzed_hits_timeseries'].empty:
                    st.success("Time-series analysis complete!")
                    st.subheader("Time-series Expression Plot for Significant Hits")
                    
                    fig = px.line(
                        st.session_state['rnaseq_analyzed_hits_timeseries'],
                        x=rnaseq_time_col,
                        y=rnaseq_value_col,
                        color='official_symbol',
                        markers=True,
                        title=f"Expression of Significant CRISPR Hits Over Time"
                    )
                    st.plotly_chart(fig, use_container_width=True)
                    st.subheader("Filtered RNA-seq Data")
                    st.dataframe(st.session_state['rnaseq_analyzed_hits_timeseries'])
                else:
                    st.warning("No RNA-seq data found for the significant CRISPR hits.")
    
    elif analysis_type == 'Polarization State Comparison':
        rnaseq_gene_col = st.selectbox("Gene Name Column:", rnaseq_cols, index=0, key="rnaseq_gene_col_pol")
        
        if st.button("Run Polarization Analysis", key="run_polarization_btn"):
            with st.spinner("Analyzing RNA-seq polarization data..."):
                st.session_state['rnaseq_analyzed_hits_polarization'] = data_analysis.analyze_rnaseq_polarization(
                    st.session_state['rna_seq_data'],
                    rnaseq_gene_col,
                    st.session_state['gene_list_for_enrichment']
                )

                if not st.session_state['rnaseq_analyzed_hits_polarization'].empty:
                    st.success("Polarization analysis complete!")
                    st.subheader("Expression Comparison by Polarization State for Significant Hits")
                    
                    fig = px.box(
                        st.session_state['rnaseq_analyzed_hits_polarization'],
                        x='Polarization_State',
                        y='Expression_Value',
                        color='Polarization_State',
                        title=f"Expression of Significant CRISPR Hits by Polarization State"
                    )
                    st.plotly_chart(fig, use_container_width=True)
                    st.subheader("Filtered and Reshaped RNA-seq Data")
                    st.dataframe(st.session_state['rnaseq_analyzed_hits_polarization'])
                else:
                    st.warning("No RNA-seq data found for the significant CRISPR hits.")

else:
    st.info("Please upload your CRISPR screen data to enable this section.")