# ETfinder

ETfinder is a framework for discovering and prioritizing RecET recombineering systems in non-model microbial chassis.  
It integrates SSB C-terminal–guided compatibility filtering, phylogeny-aware candidate ranking, genome-context analysis, and an offline web interface for interactive visualization.

---

# 1. Installation

## 1.1 Clone the repository

```bash
git clone https://github.com/lvdongyuan/ETfinder.git
cd ETfinder
```

---

# 2. Python environment

ETfinder requires **Python 3.8–3.10**  
(ETE3 is not compatible with Python ≥3.11.)

We strongly recommend installing ETfinder inside an isolated **conda environment**:

```bash
conda create -n etfinder python=3.10
conda activate etfinder
```

Then install Python dependencies:

```bash
pip install -r requirements.txt
```

`requirements.txt` contains:

```
ete3>=3.1.0
flask>=2.3.0
```

---

# 3. Install external tools (required)

ETfinder relies on the following command-line tools:

- CD-HIT  
- MAFFT  
- IQ-TREE 2/3  

Install via conda:

```bash
conda install -c bioconda cd-hit mafft iqtree
```

Verify installation:

```bash
cd-hit --version
mafft --version
iqtree --version
```

---

# 4. Initialize ETE3 taxonomy database (required once)

ETE3 requires a local NCBI taxonomy SQLite database.  
Run:

```bash
python -c "from ete3 import NCBITaxa; NCBITaxa().update_taxonomy_database()"
```

## Offline installation (no VPN access)

Place a pre-built `taxa.sqlite` file into:

```
~/.etetoolkit/taxa.sqlite
```

ETE3 will automatically detect and use this offline database.

---

# 5. Running ETfinder2

## 5.1 Offline web interface

```bash
python server/server.py
```

Open in your browser:

```
http://127.0.0.1:5050
```

Features:

- SSB-guided RecT candidate search  
- Interactive RecT & host phylogenetic trees  
- Genomic context viewer  
- Export of CSV / JSON / Newick trees  

Runs fully offline.

---

## 5.2 Running modules individually (optional)

```bash
python M1_ssb_search/ssbc7_search.py
python M2_context_tree/M2_taxontree.py
python M3_scoring/recT_scoring_new.py
```

---

# 6. Output files

After running ETfinder, results will be stored in the output directory.

### **Recommended final result**

```
Results.ctx.tsv
```

This file contains the **final ranked list of RecT candidates** after scoring and stratified prioritization.  
You may simply download this file from the web interface or locate it in the `output/` directory to inspect ranked candidates directly.

---

# 7. RecTdb: RecT–SSB Database

RecTdb includes:

- >25,000 RecT homologs  
- Associated SSB sequences  
- SSB-C-terminal motifs  
- Host taxonomy  
- Optional genomic context annotations  

---

# 8. Citation

If you use ETfinder2 or RecTdb, please cite:

> [Manuscript currently under review. Citation will be updated upon publication.]


---

# 9. License

```
The MIT License (MIT)
Copyright (c) 2025 DongyuanLv
```

---

# 10. Contact

```
DongyuanLv
ECUST
Y10240045@mail.ecust.edu.cn
```
