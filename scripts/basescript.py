import os
from pyopenms import *


# Set searchfile and database
searchfile = "/vol/VMvol/data/raw_data/peptideatlas/ftp.peptideatlas.org/distro/mzML/JD_06232014_sample1-A.mzML"
database = "/vol/VMvol/data/raw_data/peptideatlas/ftp.peptideatlas.org/distro/fasta/iPRG2015.TargDecoy.fasta"

# Set path to Percolator and PercolatorAdapter
perc_path = "/home/ubuntu/miniconda3/envs/openmsenv/bin/percolator"
percadapter_path = "/home/ubuntu/miniconda3/envs/openmsenv/bin/PercolatorAdapter"


# Run SimpleSearchEngineAlgorithm, store protein and peptide idsls

protein_ids = []
peptide_ids = []

# Enable additional score annotations
simplesearch = SimpleSearchEngineAlgorithm()
params = simplesearch.getDefaults()

score_annot = [b'fragment_mz_error_median_ppm', b'precursor_mz_error_ppm']
params.setValue(b'annotate:PSM', score_annot)

#print(params.items())
simplesearch.setParameters(params)

simplesearch.search(searchfile, database, protein_ids, peptide_ids)

# Save to disc in order to run Percolator on it
IdXMLFile().store("SSE_results.idXML", protein_ids, peptide_ids)

# Run PercolatorAdapter as command using os.system()
#print("\nStart PercolatorAdapter ... \n")
os.system(percadapter_path + " -in SSE_results.idXML -out SSE_results_percolated.idXML -percolator_executable "
          + perc_path + " -out_pin SSE_results_percolated.tab -weights SSE_results_percolated.weights "
                        "-train_best_positive -score_type q-value -generic_feature_set")


# Load the new ids
perc_protein_ids = []
perc_peptide_ids = []
IdXMLFile().load("SSE_results_percolated.idXML", perc_protein_ids, perc_peptide_ids)


# Count #hits
c_hits = 0
for p in perc_peptide_ids:
    for h in p.getHits():
        c_hits += 1


# Load fasta file
fasta_file = FASTAFile()
fasta_entries = []
fasta_file.load(database, fasta_entries)

#print("\nStart Indexing ... \n")

# Annotate the peptides with proteins
# Step is omitted, because no param "tolerance" is given
#indexer = PeptideIndexing()
#indexer_params = indexer.getParameters()
#indexer_params.setValue("tolerance", 0.05)
#indexer.setParameters(indexer_params)
#indexer.run(fasta_entries, perc_protein_ids, perc_peptide_ids)


# Annotate q-value
fdr = FalseDiscoveryRate()
fdr.apply(perc_peptide_ids)

# Filter by 1% PSM FDR
idfilter = IDFilter()
idfilter.filterHitsByScore(perc_peptide_ids, 0.01)
idfilter.removeDecoyHits(perc_peptide_ids)


# Store filtered ids
IdXMLFile().store("result_idx_fdr001.idXML", perc_protein_ids, perc_peptide_ids)


# Count #hits after filtering
c_filtered = 0
for p in perc_peptide_ids:
    for h in p.getHits():
        c_filtered += 1

print("\n")
print("Hits before filtering:", c_hits)
print("Hits after filtering:", c_filtered)


