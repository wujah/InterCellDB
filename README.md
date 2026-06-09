
This repository contains two complementary Python-based toolkits designed for Google Colab and Jupyter Notebook environments to support computational analysis of protein structure interactions and cell-type interaction networks. The Protein Structure Distance and Interface Analysis Toolkit uses Biopython and Pandas to parse PDB files, calculate atom-, residue-, and chain-level distance metrics, and evaluate overlap between protein interaction interfaces. The Cell Type Interaction Network & Community Detection Toolkit processes pairwise Jaccard similarity data to construct weighted interaction graphs, identify communities using a customized weighted Clauset-Newman-Moore algorithm, and visualize highly similar cell-type interaction networks. Together, these tools provide a flexible workflow for investigating molecular interfaces, interaction overlap, and higher-order network organization across biological systems.


<img width="335" height="438" alt="toolkit1" src="https://github.com/user-attachments/assets/1d6b0a9a-6e8a-4317-a421-d413262c8073" />


Protein Structure Distance and Interface Analysis Toolkit

Description
This Python program, designed for Google Colab or Jupyter Notebook environments, provides a suite of tools for analyzing protein structures. 
It utilizes Biopython to parse Protein Data Bank (PDB) files and calculates various distance metrics between chains, residues, and atoms to identify potential interfaces and interactions. 
Additionally, it processes external interaction data using Pandas to determine interface overlaps between different protein chains or components.


Prerequisites:
Ensure the following libraries are installed in your Python environment:
biopython (Handles PDB parsing and spatial data extraction)
pandas (Handles flat-file data processing)
math (Standard library for coordinate mathematics)


Core Features:
PDB Parsing: Uploads and directly parses structural data, extracting coordinates for specific models, chains, residues, and atoms.
Inter-Atomic Distances: Calculates the Euclidean distance between individual atoms across different protein chains.
Residue-Level Distances: Computes minimal distances between all atoms of specific residues to identify close contacts (< 3.0 Å).
Cα Backbone Mapping: Extracts C-alpha (Cα) coordinates to compute backbone-level distances between entire chains.
Interface Overlap Analysis: Reads flat .dat files detailing interacting pairs and computes an overlap matrix to determine if specific components share common interaction interfaces.


Usage Instructions:
Step 1: Environment Setup
Run the first block to execute !pip install biopython to ensure structural parsing dependencies are met.
Step 2: Upload PDB File
Execute the file upload cell. When prompted by the Colab widget, upload your target .pdb file. The program will initialize the parser and output a hierarchy of the models, chains, and atomic coordinates.
Step 3: Calculate Distances
Run the respective distance calculation cells depending on your analytical needs. You can manually adjust the chain_id_protein1 and chain_id_protein2 variables in the code blocks to target specific chains 
(e.g., 'A' and 'B').
Step 4: Upload Interaction Data
Execute the second file upload cell to input your interface data (e.g., FGFR1_InterfaceOverlap.dat).
Step 5: Analyze Interface Overlap
Run the final Pandas and looping cells. The code utilizes delim_whitespace=True to parse the flat file regardless of spacing inconsistencies. It will cross-reference the interaction nodes (columns 'I' and 'J') 
and print a binary matrix indicating overlap status (1 for overlap, 0 for no overlap).


Notes and TroubleshootingData Formatting: 
If the Pandas parser encounters errors with your .dat file, verify that the file contains consistent columns without trailing unstructured text.
Heavy Computations:
Calculating the distance from every atom in Chain A to every atom in Chain B has a time complexity of O(N*M). 
For very large macro-complexes, extracting Cα distances (Cell 6) or specific residue distances (Cell 5) is significantly faster than the all-atom approach.


<img width="356" height="391" alt="toolkit2" src="https://github.com/user-attachments/assets/876ac6c3-fdc0-4d3b-8a28-f7151a4ec175" />


Cell Type Interaction Network & Community Detection Toolkit

Description:
This Python program is designed for Google Colab or Jupyter Notebook environments to analyze and visualize the similarity of interactions between different cell types. It processes pairwise Jaccard similarity data, constructs a weighted network graph, and identifies structural communities using a customized Clauset-Newman-Moore (CNM) agglomerative algorithm adapted for weighted edges.

Prerequisites:
Ensure the following Python libraries are installed in your environment:
pandas (For data manipulation and flat-file parsing)
networkx (For graph construction and network operations)
matplotlib (For network visualization)
numpy (For numerical operations and infinity representations)
google.colab (If running in Colab, for the file upload widget)

Data Requirements
The script expects a flat data file named PPICellAdhesion_InterCellType_JaccardSimilarity.dat.
The file should contain three whitespace-separated columns:
I: The integer ID of the first cell type.
J: The integer ID of the second cell type.
Similarity: A float representing the Jaccard similarity score between cell type I and cell type J.


Core Workflow:
1. Data Ingestion
Uses Colab's files.upload() widget to import the .dat file.
Parses the file using pandas.read_csv() with flexible whitespace delimiting to gracefully handle irregular spacing.

2. Graph Construction
Initializes a networkx Graph object.
Iterates through the data to establish nodes and edges.
Filtering Threshold: By default, only cell type pairs with a Similarity > 0.25 are added as 
connected edges in the graph, with the similarity score serving as the edge weight.

3. Community Detection (Weighted CNM Algorithm)
The script implements a custom version of the Clauset-Newman-Moore algorithm:
Initialization: Every node starts in its own isolated community.
Modularity Optimization ($Q$): Evaluates potential community merges by calculating the change in network modularity. The custom modularity_weighted function accounts for the continuous edge weights rather than just binary connections.
Agglomeration: Iteratively merges the pair of communities that yields the highest increase in modularity until no further improvement can be made.

4. Visualization
Generates a 2D network plot using matplotlib and networkx.draw().
Nodes are represented as sky-blue circles with bold labels, visually mapping the highly similar (>0.25) cell interactions.


Usage Instructions
Run the initial cell to trigger the upload widget and select your PPICellAdhesion_InterCellType_JaccardSimilarity.dat file.
Execute the Pandas parsing cell to format the arrays.
Run the main executable block to build the network, compute the CNM communities, and print the community cluster assignments to the console.
Run the final block to render the topological visualization of the network.

Customization
Threshold Adjustment: 
To analyze weaker or stronger interaction networks, modify the if similarity > 0.25: line in the graph generation loops.
Modularity Tuning:
The alpha parameter in the modularity_weighted function (currently alpha = 0.5) can be adjusted to penalize or favor intra-community edge weights versus expected degree sums differently.
