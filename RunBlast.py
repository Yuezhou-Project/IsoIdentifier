
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from os import path


def globalRun (d_dataset, p_dir_blast, debug = 1) : 
    
    for PDB_ID in d_dataset.keys () : 
        if d_dataset[PDB_ID]["conserve"] == 1 : 
            p_fasta = d_dataset[PDB_ID]["best"]["fasta"]
            p_out_blast = p_dir_blast + PDB_ID + ".xml"
            blastp_path = "/home/buhan/Desktop/soft/ncbi-blast-2.10.1+/bin/blastp"#change
            blastp_cline = NcbiblastpCommandline(cmd=blastp_path,query=p_fasta, db="/home/buhan/Desktop/myproject/pdb", outfmt=5, out=p_out_blast)#change
            if debug : print(blastp_cline)
            if not path.exists(p_out_blast) : 
                stdout, stderr = blastp_cline()
            d_dataset[PDB_ID]["xml"] = p_out_blast
            d_dataset[PDB_ID]["align"] = {}
            # parse blast out
            result_handle = open(p_out_blast)
            blast_records = NCBIXML.read(result_handle)
            for alignment in blast_records.alignments:
                for hsp in alignment.hsps:
#                     print alignment.title
                    PDB_find = alignment.title.split ("|")[4].split (" ")[0]
                    d_dataset[PDB_ID]["align"][PDB_find] = hsp.expect
            
            result_handle.close ()





