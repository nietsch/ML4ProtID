from pyopenms import *

# Set searchfile and database
searchfile = "/vol/VMvol/data/raw_data/peptideatlas/ftp.peptideatlas.org/distro/mzML/JD_06232014_sample1-A.mzML"
database = "/vol/VMvol/data/raw_data/peptideatlas/ftp.peptideatlas.org/distro/fasta/iPRG2015.TargDecoy.fasta"

# Set path to Percolator and PercolatorAdapter
perc_path = "/home/ubuntu/miniconda3/envs/openmsenv/bin/percolator"
percadapter_path = "/home/ubuntu/miniconda3/envs/openmsenv/bin/PercolatorAdapter"


# Run SimpleSearchAlgorithm, store protein and peptide ids
protein_ids = []
peptide_ids = []

simplesearch = SimpleSearchEngineAlgorithm()
params = simplesearch.getDefaults()

# Enable additional scroe annotations
score_annot = [b'fragment_mz_error_median_ppm', b'precursor_mz_error_ppm']
params.setValue(b'annotate:PSM', score_annot)

params.setValue(b'annotate:PSM', score_annot)

simplesearch.setParameters(params)

simplesearch.search(searchfile, database, protein_ids, peptide_ids)


# Count #hits
c_hits = 0
for p in peptide_ids:
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


# Without Percolator -> less hits
# Annotate q-value
fdr = FalseDiscoveryRate()
fdr.apply(peptide_ids)

# Filter by 1% PSM FDR
idfilter = IDFilter()
idfilter.filterHitsByScore(peptide_ids, 0.01)
idfilter.removeDecoyHits(peptide_ids)

c = 0
for p in peptide_ids:
    for h in p.getHits():
        c = c + 1


# Store filtered ids
IdXMLFile().store("noperc_result_idx_fdr001.idXML", protein_ids, peptide_ids)

print("\n")
print("Hits before filtering:", c_hits)
print("Hits after filtering (without Percolator):", c)




