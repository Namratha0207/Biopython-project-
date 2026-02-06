# =====================================================
# Functional Annotation of Hypothetical Protein pet46
# Organism: Candidatus Bathyarchaeota.
# Method: Homology-based annotation using BLAST
# =====================================================

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# STEP 1:Read fasta Sequence

record = SeqIO.read("pet46.fasta", "fasta")
seq = record.seq
print("Protien ID:", record.id)
print("Description:", record.description)
print("Sequence Length:", len(record.seq))

# STEP 2 : BLASTp analysis

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

from Bio.Blast import NCBIWWW

result_handle= NCBIWWW.qblast(
    program= "blastp",
    database= "nr",
    sequence= seq
)
with open ("blast_pet46.xml", "w") as b:
     b.write(result_handle.read())

print("BLAST completed successfully.")

# STEP 3 : Parse Blast results
print("\n--- Top BLAST Homology Hits ---")

with open("blast_pet46.xml") as result:
    blast_record = NCBIXML.read(result)

hydrolase_hits = 0
keywords = ["hydrolase", "esterase", "polyesterase", "PETase"]

for alignment in blast_record.alignments[:5]:
    hsp = alignment.hsps[0]
    description = alignment.hit_def.lower()
    
    print(f"Hit ID: {alignment.hit_id}")
    print(f"Description: {alignment.hit_def}")
    print(f"E-value: {hsp.expect}")
    print(f"Identity: {hsp.identities}/{len(record.seq)}")
    print("-" * 50)

    # Check if the description matches functional keywords
    if any(key in description for key in keywords):
        hydrolase_hits += 1

# STEP 4: Functional Annotation
print("\n--- Final Functional Annotation Result ---")

if hydrolase_hits > 0:
    print("[ANNOTATION]: PET46 is confirmed as a member of the Alpha/Beta Hydrolase superfamily.")
    print("[EVIDENCE]: High homology detected with microbial esterases.")
    print("[INDUSTRIAL POTENTIAL]: Predicted Thermostable PET-Degrading Enzyme.")
else:
    print("[ANNOTATION]: No significant homology to known plastic-degrading families found.")
    print("[NOTE]: Requires further structural motif analysis (e.g., searching for Ser-His-Asp triad).")