#Codes to identify local structural replacement of glycosyl moiety of ligands from open-accessible data repository Protein Data Bank
### Dependencies  ###
###################

TM align -> path l.13 script RunOtherSoft.py // http://zhanglab.ccmb.med.umich.edu/TM-align/

ShaEP -> path l.14 script RunOtherSoft.py // http://users.abo.fi/mivainio/shaep/

blastp

Requirements (Centos 7): 
- openbabel
- EMBOSS
- R 

Python (3.6) modules:
- os
- copy
- re
- subprocess
- numpy
- biopython (Bio.blast)
- shutil
- urllib
- math



### Scripts py descriptions  ###
###########################

- RunBlast.py: run Blast queries
- analysis.py: automated statistical analysis and draw R plots
- arrangeResult.py: manage folders containing LSR
- buildData.py: extract files from the PDB
- cleanResult.py: clean folder result, remove temp files
- downloadFile.py: Import fasta file and pdb file using urlretrieve protocol
- ionSearch.py: search and analysis metals in the binding site
- managePDB.py: clean and decompress PDB extracted
- neighborSearch.py: define binding site
- parseEMBOSS.py: parsing EMBOSS output
- parsePDB.py: parsing pdb file
- parseShaep.py: parsing ShaEP output
- parseTMalign.py: parsing TMalign ouput
- pathManage.py: path manager
- refClassification.py: build pdb classification
- runOtherSoft.py: run external tools
- smileAnalysis.py: assign SMILES categories
- substructTools.py: build reference substructure
- superimpose.py: superimpose ligand using Kabsh's algorithm
- superposeStructure.py: superimpose protein
- tool.py: tool box functions
- writePDBfile.py: write pdb files

- main.py: MAIN with parameters

thresold_RX = 2.7: minimal structure resolution
thresold_BS = 4.5: distance threshold to define the binding site from ligand
thresold_blast = 1e-100: threshold blast
thresold_superimposed_ribose = 2.5: distance threshold to extract LSR from ligand
thresold_IDseq = 100: threshold identity sequence for similar PDB
thresold_shaep = 0.2: threshold of shaep overlap
l_ligand_out = ["AMP", "ADP", "ATP", "ACP", "AD9", "NAD", "AGS", "U5P", "UDP","UTP", "APC", "C5P","CDP","CTP", "AOV", "ANP","5GP", "GDP", "GTP", "ANP"]: list of ligand not considered

### other scripts ###
##################
-downloadPDB.sh: download PDB database

# R scripts for plots
- histogram.R
- histograms.R
- histogramsRMSD.R
- barplotQuantity.R
- piePlot.R
- piePlot_count.R

### Paths management ###
########################

BLAST in local using PDB sequences ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz

# 1. create db for blast #
#####################

# format ncbi format
seqret pdb.fasta pdb_ncbi.fasta -osf ncbi

# make db
makeblastdb -in pdb_ncbi.fasta -out pdb -dbtype prot

# 2. run query blast for test fasta
###################################

blastp -query test.fasta -db pdb -out result.txt
# 3. run main.py
########################

This is only code need to be run, then it will invoke other codes to automatly parse the PDB.

# 4. Previous codes are referred to construct current ones.

https://github.com/ABorrel/LSRs
