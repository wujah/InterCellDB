# InterCellDB

This repository contains two complementary Python notebooks for analyzing protein-binding interfaces and similarities among cell–cell adhesion interaction networks:

- **CalculateProteinInterfaceNew.ipynb** calculates distances between protein chains, identifies contact residues, exports interface residues, and compares interface overlap among binding partners.
- **InterCellDB.ipynb** reads precomputed pairwise Jaccard similarity scores, constructs a weighted network of cell–cell interfaces, identifies communities using a custom greedy merging routine, and visualizes the network.

Both notebooks use Google Colab file-upload widgets. Google Colab is the simplest way to run them. To use a local Jupyter environment, replace the Colab upload/download calls with local file paths and file handling.

## Contents

- [Protein structure distance and interface analysis](#protein-structure-distance-and-interface-analysis)
- [Cell–cell interface network and community detection](#cellcell-interface-network-and-community-detection)

## Protein structure distance and interface analysis

**Notebook:** `CalculateProteinInterfaceNew.ipynb`

<img width="335" height="438" alt="Protein structure distance and interface analysis workflow" src="https://github.com/user-attachments/assets/1d6b0a9a-6e8a-4317-a421-d413262c8073" />

### Description

This notebook uses Biopython to read a protein-complex PDB file and calculate distances between atoms and residues in two selected chains. It identifies contacting residue pairs using a minimum interatomic distance of less than 3.0 Å and includes a cell that exports interface residues to CSV. An optional workflow compares the target-protein residues contacted by different binding partners to identify interface overlap.

### Prerequisites

- **Biopython:** PDB parsing and spatial searches.
- **pandas:** Reading the optional interface-overlap table.
- **Python standard library:** `math`, `csv`, and `pathlib`; no separate installation is needed.
- **Google Colab:** File-upload and download widgets.

In Colab, run the first cell to install Biopython. Install pandas if it is not already available:

```python
!pip install biopython pandas
```

### Input files

#### Protein-complex structure

Upload one `.pdb` file containing both interacting proteins as separate chains in the same coordinate system. The default chain identifiers are `A` and `B`. Separate, unaligned monomer structures are not sufficient for interface analysis.

Use a single-model structure for the original distance-analysis cells, which do not separate coordinates by model. The CSV export cell explicitly uses the first model.

#### Optional interface-overlap table

The overlap workflow expects the filename `FGFR1_InterfaceOverlap (3).dat`. If your file has another name, update the filename in `pandas.read_csv()`.

The file must be whitespace-delimited and contain these headers:

| Column | Meaning |
| --- | --- |
| `I` | Integer identifier of a binding partner. |
| `J` | Target-protein residue number contacted by that partner. |

For example:

```text
I J
1 25
1 26
2 26
2 80
```

Partners 1 and 2 share target residue 26 in this example. All partner records must refer to the same target protein and use consistent residue numbering. Prepare this table separately; the exported interface-residue CSV is not automatically converted into this format.

### Usage instructions

Run the required cells in notebook order. Cell numbers below count all notebook cells, including the Markdown cell above the CSV export.

#### Step 1 — Set up the environment

Run Cell 1 to install Biopython.

#### Step 2 — Upload and inspect the PDB file

Run Cell 2 and select one complex `.pdb` file. Inspect the printed model and chain information to confirm the chain identifiers.

#### Step 3 — Calculate distances and identify contacts

Run the relevant cells for your analysis:

| Cell | Purpose |
| --- | --- |
| 3 | Demonstrate distances from one selected atom to atoms in the other chain. |
| 4 | Calculate all cross-chain atom-pair distances. |
| 5 | Calculate atom-pair distances between two specified residues. |
| 6 | Calculate cross-chain Cα distances. |
| 7 | Calculate the minimum atom-pair distance for every cross-chain residue pair. |
| 8 | Print residue pairs with minimum atom distance strictly less than 3.0 Å. |

Update `chain_id_protein1` and `chain_id_protein2` in each cell you run if your chains are not `A` and `B`. For Cell 5, also update the residue numbers to match the PDB numbering. Run a cell importing `math` before any distance cell that uses it.

#### Step 4 — Export interface residues

Run Cell 10, headed **Save and download the interface-residue CSV**, after loading the structure and setting the chain identifiers. The cell independently calculates contacts and downloads `<PDB_filename_stem>_interface_residues.csv`.

The default cutoff is strictly less than 3.0 Å. To change the export cutoff, edit `cutoff_angstrom` in this cell. Changing the cutoff in Cell 8 does not change the export cutoff.

#### Step 5 — Upload the optional overlap table

To compare different partners, run Cell 11 and upload the interface-overlap `.dat` file. Skip Cells 11–14 if you only need contacts from one complex.

#### Step 6 — Compare partner interfaces

Run Cells 12 and 13 to read the table and select columns `I` and `J`. Then run Cell 14 to print the overlap status for each partner pair. The current loops compare partner IDs 1–8 and target residue numbers 0–349; adjust the ranges for other datasets.

### Outputs

- Printed atom- and residue-distance summaries, in Å.
- Printed contacting residue pairs satisfying the selected distance criterion.
- A downloaded CSV containing one row per interface residue from either chain, with chain identifier, residue number, insertion code, residue name, hetero flag, and minimum distance to the other chain.
- Optional printed partner-pair overlap flags: `1` indicates a shared target residue within the tested range, and `0` indicates none was found in that range.

### Notes and troubleshooting

#### Data formatting

Check that `.dat` files contain the required headers and a consistent number of whitespace-separated fields. If parsing fails, correct the input before continuing. The overlap-table parsing cell contains an existing diagnostic that refers to the similarity-table filename; that filename is not an additional input requirement for the overlap workflow.

#### Computation time

Comparing every atom in one chain with every atom in another requires a number of distance evaluations proportional to the product of their atom counts. For large complexes, selected-residue or Cα calculations are less expensive, but Cα distances are not equivalent to minimum all-atom contact distances.

#### Contact interpretation

The distance calculations include all atoms present in the selected chains, including hydrogen atoms and non-protein residues when supplied. The original residue-distance cells do not distinguish insertion codes in their residue-number lookups. Use unambiguous residue numbering for these cells.

Shared contact residues suggest interface overlap but do not prove competitive binding. Conversely, non-overlapping residue lists do not establish that two partners can bind simultaneously.

## Cell–cell interface network and community detection

**Notebook:** `InterCellDB.ipynb`

<img width="356" height="391" alt="Cell–cell interface similarity network workflow" src="https://github.com/user-attachments/assets/876ac6c3-fdc0-4d3b-8a28-f7151a4ec175" />

### Description

This notebook analyzes a supplied table of pairwise Jaccard similarity scores between cell–cell adhesion interaction fingerprints. Each network node represents a cell–cell interface, defined by a pair of interacting cell types. Edges connect interfaces whose similarity is greater than 0.25, with the similarity score used as the edge weight.

A custom greedy community-detection routine evaluates possible community merges using the implemented `modularity_weighted()` score and continues while that score improves. The final cell draws the thresholded similarity network. The notebook does not construct fingerprints or calculate their Jaccard similarity scores.

### Prerequisites

- **pandas:** Reading and inspecting the similarity table.
- **NetworkX:** Network construction and drawing.
- **Matplotlib:** Plot display.
- **NumPy:** Numerical utilities.
- **Google Colab:** File-upload widget.

Install missing packages before running the notebook:

```python
!pip install pandas networkx matplotlib numpy
```

### Input files

Upload `PPICellAdhesion_InterCellType_JaccardSimilarity.dat`. If your filename differs, update it in `pandas.read_csv()`.

The file must contain three whitespace-separated columns:

| Column | Meaning |
| --- | --- |
| `I` | Integer identifier of the first cell–cell interface. |
| `J` | Integer identifier of the second cell–cell interface. |
| `Similarity` | Pairwise Jaccard similarity score, ranging from 0 to 1. |

For example:

```text
I J Similarity
1 2 0.30
1 3 0.12
2 3 0.40
```

This example illustrates the format only. The current notebook loops inspect 820 rows and scan `I` from 0 to 40 and `J` from `I` to 41. Review these limits and your identifier convention before using a different dataset; a three-row example cannot be run unchanged through these loops.

### Usage instructions

Run the cells in notebook order.

#### Step 1 — Upload the similarity table

Run Cell 1 and select the `.dat` file.

#### Step 2 — Read and inspect the data

Run Cell 2 to create `dfX` and inspect the table preview and column information. Confirm the exact headers `I`, `J`, and `Similarity`.

#### Step 3 — Select the analysis columns

Run Cell 3 to create and print the interface identifier and similarity columns.

#### Step 4 — Inspect the interface pairs

Cell 4 prints the pairs and scores matching the hard-coded ranges. This optional inspection cell has an existing leading-indentation error; remove the extra indentation before the outer loop and its body before running it. Later cells do not depend on its output.

#### Step 5 — Build the network and detect communities

Run Cell 5 to define the scoring and merging functions, construct the network, and print community memberships. The implementation is custom and should not be assumed equivalent to NetworkX's standard CNM implementation.

#### Step 6 — Visualize the network

Run Cell 6 to rebuild and draw the thresholded network. Node labels are interface identifiers. All nodes use the same color; the plot does not encode community assignments.

### Outputs

- A preview of the input similarity table and selected columns.
- Printed community memberships for the constructed network.
- An on-screen network plot of interface pairs with similarity greater than 0.25.

The notebook does not automatically save community assignments or the plot to an output file.

### Notes and troubleshooting

#### Data formatting

Use numeric interface identifiers and similarity scores. Correct any parsing errors before running subsequent cells. The parsing diagnostic prints line 227 when a parser error occurs; it does not recover a valid table automatically.

#### Dataset size and identifier ranges

The loops use hard-coded limits for the original dataset. Update the row and identifier ranges in both network-construction cells when analyzing a different table. Interfaces without retained edges are not explicitly added, so they do not appear as isolated nodes or communities.

#### Similarity threshold

To change the edge cutoff, edit `if similarity > 0.25:` in both Cells 5 and 6. Keep the values consistent so the community analysis and visualization use the same network.

#### Community scoring

The custom scoring function sets `alpha = 0.5`. Changing this value changes the implemented scoring criterion and may change the detected communities. Interpret results according to the formula in the notebook.

#### Plot layout

The drawing does not fix a layout seed, so node positions can differ between runs even when the graph is unchanged.
