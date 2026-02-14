# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    aligner = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1) # create aligner
    
    # store species names and alignment scores
    species_scores = []
    
    # align each species to human
    gg_score, _, _ = aligner.align(hs_seq, gg_seq)
    species_scores.append(("Gallus gallus", gg_score))
    
    mm_score, _, _ = aligner.align(hs_seq, mm_seq)
    species_scores.append(("Mus musculus", mm_score))
    
    br_score, _, _ = aligner.align(hs_seq, br_seq)
    species_scores.append(("Balaeniceps rex", br_score))
    
    tt_score, _, _ = aligner.align(hs_seq, tt_seq)
    species_scores.append(("Tursiops truncatus", tt_score))
    
    # in order of most similar to human BRD - as instructed above
    species_scores.sort(key=lambda x: x[1], reverse=True)
    
    # print results in order of similarity
    print("Species ordered by similarity to Homo sapiens BRD2:")
    for species, score in species_scores:
        print(f"{species}: {score}")

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    print("\nAlignment scores for each species with Homo sapiens BRD2:")
    print(f"Gallus gallus: {gg_score}")
    print(f"Mus musculus: {mm_score}")
    print(f"Balaeniceps rex: {br_score}")
    print(f"Tursiops truncatus: {tt_score}")
    

if __name__ == "__main__":
    main()
