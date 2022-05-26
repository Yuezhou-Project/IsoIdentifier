"""
BORREL Alexandre
04-2012
"""
# global module
from urllib.request import urlretrieve
from os import system, path, listdir
from re import search
from shutil import copy

# personal module
     

def importPDB ( PDB_ID , directory, dir_by_PDB = 1, debug=0 , dbPDB = "" ):
    """Retrieve in http://www.pdb.org PDB file
    args: - list_PDB
          - directory out
    return: - error
            - write in log directory errors download
    NB: make log file best
    """
    name_PDB = PDB_ID.upper () + ".pdb"
    if dir_by_PDB : 
        pr_out = directory + PDB_ID + "/"
        if not path.isdir(pr_out) : 
            system ( "mkdir " + pr_out )
    else : 
        pr_out = directory
    
    if dbPDB != "" : 
        
        if path.exists(pr_out + PDB_ID.upper () + ".pdb") and path.getsize(pr_out + PDB_ID.upper () + ".pdb") != 0: 
            return pr_out + PDB_ID.upper () + ".pdb"
        else : 
            if not  path.exists(dbPDB + PDB_ID.lower() + ".pdb") : 
                return 0
            else : 
                copy (dbPDB + PDB_ID.lower() + ".pdb", pr_out + PDB_ID.upper () + ".pdb")
                return pr_out + PDB_ID.upper () + ".pdb"
        return 0

    else : 
    #     print PDB_ID, "Download"
        adresseSeq = ( "http://www.pdb.org/pdb/files/%s.pdb" % PDB_ID )
        try:
            name_PDB = PDB_ID.upper () + ".pdb"
            if dir_by_PDB : 
                pr_out = directory + PDB_ID + "/"
                if not path.isdir(pr_out) : 
                    system ( "mkdir " + pr_out )
            else : 
                pr_out = directory
            
            p_filout = pr_out + name_PDB
            if path.exists(p_filout) : 
                return p_filout
            else : 
                path_file_pdb = urlretrieve( adresseSeq )
                if debug : print(path_file_pdb)
                cmd = "mv " + path_file_pdb[0] + " " + pr_out + name_PDB
                if debug : print(cmd)
                system ( cmd )
                print(str( PDB_ID ) + "-> done")
                return pr_out + name_PDB
        except:
            print(str( PDB_ID ) + "-> ERROR DOWNLOAD PDB file")
            return 0
    
    
def importFasta ( PDB_ID ,directory, dir_by_PDB = 1, debug=0, fastaGlobal = "" ):
    """Retrieve in http://www.pdb.org PDB file
    args: - list_PDB
          - repertory out
    return: - NULL (print error)
    """
    name_PDB = PDB_ID.upper () + ".fasta"
    if dir_by_PDB : 
        pr_out = directory + PDB_ID + "/"
        if not path.isdir(pr_out) : 
            system ( "mkdir " + pr_out )
    else : 
        pr_out = directory
        
    if fastaGlobal != "" : 
        
        if path.exists(pr_out + name_PDB) and path.getsize(pr_out + name_PDB) != 0: 
            return pr_out + name_PDB
        else :
            filout = open (pr_out + name_PDB, "w")
            filin = open (fastaGlobal, "r")
            l_lines_fasta = filin.readlines()
            filin.close ()
            nb_lines = len (l_lines_fasta)
            i = 0
            while i < nb_lines : 
                if len (PDB_ID) == 4 : 
                    if search (">" + PDB_ID.lower(), l_lines_fasta[i]) : 
                        filout.write (l_lines_fasta[i])
                        filout.write (l_lines_fasta[i + 1])
                        i = i + 2
                    else : 
                        i = i + 1
                        
                else : 
                    if search (">" + PDB_ID, l_lines_fasta[i]) : 
                        filout.write (l_lines_fasta[i])
                        filout.write (l_lines_fasta[i + 1])
                        i = i + 2
                    else : 
                        i = i + 1
                    
            filout.close ()
            return pr_out + name_PDB
            
    else : 
    #     print PDB_ID, "download Fasta"
        adresseSeq = ( "http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=%s.pdb" % PDB_ID )
        
        try:
            p_filout = pr_out + name_PDB
            if path.exists(p_filout) :    
                return p_filout
            else : 
                path_file_pdb = urlretrieve( adresseSeq )
                if debug : print(path_file_pdb)
                cmd = "mv " + path_file_pdb[0] + " " + pr_out + name_PDB
                if debug : print(cmd)
                system ( cmd )
                print (str( PDB_ID ) + "-> done")
                return pr_out + name_PDB
        except:
            print(str( PDB_ID ) + "-> ERROR DOWNLOAD FASTA file")
            return 0
