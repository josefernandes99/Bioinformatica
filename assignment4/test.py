# Import bioservices module, to run remote UniProt queries
from bioservices import UniProt
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import os

def get_protein_info(protein_id):
    service = UniProt()
    info = service.search(protein_id)
    return info

def get_protein_sequence(protein_id):
    service = UniProt()
    sequence = service.retrieve(protein_id, frmt="fasta")
    return sequence

def get_top_matching_sequences(query_sequence, database):
    # Perform BLAST search against specified database
    result_handle = NCBIWWW.qblast("blastp", database, query_sequence)
    
    # Parse XML output
    blast_records = NCBIXML.parse(result_handle)
    blast_record = next(blast_records)
    
    # Extract top matching sequences and scores for unique species (excluding human)
    top_hits = {}
    for alignment in blast_record.alignments:
        title = alignment.title
        species = None
        for part in title.split("["):
            if "Homo sapiens" not in part:
                species = part.split("]")[0].strip()
                break
        if species is None:
            continue
        sequence = alignment.hsps[0].sbjct # extract matching sequence
        score = alignment.hsps[0].score # extract matching score
        
        # Keep only the top-scoring sequence for each unique species
        if species not in top_hits or score > top_hits[species][2]:
            top_hits[species] = (species, [], score)
        
        # Keep up to 10 top-scoring sequences for each unique species
        if len(top_hits[species][1]) < 10:
            top_hits[species][1].append(sequence)
    
    # Sort the top hits by score and return the top 10 for unique species
    top_species = []
    for species, sequences, score in sorted(top_hits.values(), key=lambda x: x[2], reverse=True):
        for sequence in sequences:
            top_species.append((species, sequence))
            if len(top_species) == 10:
                return top_species
    
    return top_species


def main():
    os.system("cls") 
    
    print("\n###################     EX1     ###################")    
    #ex1
    #protein_name, protein_id = input("Enter the protein name: ").split(" ")
    #print("Protein name: " + protein_name + "\nProtein ID: " + protein_id)
    protein_id = "P68871"
    info = get_protein_info(protein_id)
    sequence = get_protein_sequence(protein_id)
    info = info.split("\n")
    info1 = info[0].split("\t")
    info2 = info[1].split("\t")
    for x in range(len(info1)):
        print(info1[x] + ": " + info2[x])
    print("\nsequence: ")
    print(sequence)
    
    #ex2
    print("###################     EX2     ###################")
    print("############ pode demorar algum tempo #############")
    #database = input("Enter the database: ")
    database = "nr"
    top_species = get_top_matching_sequences(sequence, database)
    for i, (species, sequence) in enumerate(top_species, 1):
        print(f"Match {i} - Species: {species}\nSequence: {sequence}\n")
        
    #ex3
    print("###################     EX3     ###################")
    
    
main()
