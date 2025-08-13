# CRISPR Screen Data Analysis & AI Report Generator
This Streamlit application provides a user-friendly interface to analyze CRISPR screen data, perform gene set enrichment analysis, and then generates an AI-powered research report for using Gemini's API. 

## Features
Data Upload: Upload tab-separated (.tsv) or CSV (.csv) CRISPR screen data files.
Data Preprocessing: Cleans up initial data for processing.
Interactive Filtering: Filters significant gene hits based on user-defined FDR thresholds and minimum guides per gene.
Main Enrichment Analysis: Perform gene set enrichment analysis using gProfiler across various biological sources (e.g., GO, KEGG).
AI-Powered Main Report Generation: Generates a comprehensive research report summarizing key gene hits and enriched pathways using the Gemini AI model, with an option for custom input depending on research focus.
Targeted Pathway Analysis: Further analysis into specific pathways identified in the main enrichment, generating detailed AI reports for their biological roles and contributing genes.

## Prerequisites
Before running this application, ensure you have the following installed:

* Python 3.8
* pip (Python package installer)
* Google Cloud Project with Vertex AI API enabled.
* Service Account Key for authentication (recommended for local development).

## Installation
Enter following commands into terminal:

    bash
    python -m venv venv
    source venv/bin/activate

    pip install streamlit pandas google-cloud-aiplatform gprofiler-official 


## Google Cloud & Gemini API Setup

To use the Gemini AI features, you need to configure access to Google Cloud's Vertex AI API.

* Find your Google Cloud Project ID (e.g., my-project-12345).

* Choose a Vertex AI location (e.g., us-central1).

* Create .streamlit/secrets.toml:

* In your app's root directory, create a folder named .streamlit.

* Create secrets.toml within .streamlit: 

Ini, TOML
gcp_project_id = "your-google-cloud-project-id"
gcp_location = "us-central1"
Add .streamlit/secrets.toml to your .gitignore file.

* Uncomment vertexai.init in streamlit_app.py:

* Ensure this line is active:

vertexai.init(project=st.secrets["gcp_project_id"], location=st.secrets["gcp_location"])

## To use

1. Open terminal and run:
    streamlit run streamlit_app.py
   This should open in browser

2.  Upload Data:
    * Either CSV or TSV are accepted 

3.  Filter Genes:
    * Use the filters in the sidebar to adjust 'Hit Type', 'FDR Threshold', and 'Minimum Guides Per Gene'

4.  Run Main Enrichment Analysis
    * Click the "Run Main Enrichment Analysis" button.

5.  Generate AI Report:
    * Enter information requested in prompt boxes
    * Click "Generate AI Report with Gemini" after Main Enrichment

6.  Targeted Pathway Analysis (Optional):
    * For further analysis select pathway name in the drop box and adjust FDR
    * Can optionally generate another AI analysis specifically focused on the new inputs

## Input Data Format

The application primarily expects a TSV or CSV without headers but can attempt to convert to expected tab-seperated, no-header format. 

There should be 14 columns and the required columns (and their mapping in your `data_analysis.py`) should include:
* `id` (Gene Symbol)
* `neg|fdr` (FDR for negative selection)
* `neg|lfc` (Log2 Fold Change for negative selection)
* `pos|fdr` (FDR for positive selection)
* `pos|lfc` (Log2 Fold Change for positive selection)
* `n_guides` (Number of guides per gene)

The `data_analysis.py` module is responsible for mapping the raw input columns to these standardized names if the input file uses different nomenclature.
