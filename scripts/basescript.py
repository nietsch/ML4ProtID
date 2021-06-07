import os
from CTDopts.CTDopts import CTDModel
from CTDsupport import *
from pyopenms import *


def main():

    # register command line arguments
    model = CTDModel(
        name="ML4ProtID",
        version="1.0",
        description="Script for finding peptide candidates incorporating predicted theoretical peak intensities.",
        docurl="",
        category="",
        executableName="",
        executablePath=""
    )
    # Register input mzML
    model.add(
        "input",
        required=True,
        type="input-file",
        is_list=False,
        file_formats=["mzml"],
        description="Input file"
    )

    # Register target/decoy database fasta file for database search
    model.add(
        "database",
        required=True,
        type="input-file",
        is_list=False,
        file_formats=["fasta"],
        description="Target/Decoy database file"
    )

    # Register path to Percolator executable, needed by OpenMS tool PercolatorAdapter
    model.add(
        "percolator_path",
        required=True,
        type="string",
        is_list=False,
        description="Path to Percolator executable"
    )

    # Register path to PercolatorAdapter executable in order to run it with os.system
    model.add(
        "percolator_adapter_path",
        required=True,
        type="string",
        is_list=False,
        description="Path to PercolatorAdapter executable"
    )

    model.add(
        "output",
        required=True,
        type="output-file",
        is_list=False,
        file_formats=["idxml"],
        description="Output file"
    )

    defaults = {}
    addParamToCTDopts(defaults, model)

    # parse command line
    # if -write_ini is provided, store model in CTD file, exit with error code 0
    # if -ini is provided, load CTD file into defaults Param object and return new model with parameters set as defaults
    arg_dict, openms_params = parseCTDCommandLine(sys.argv, model, defaults)

    # Set the arguments
    searchfile = arg_dict["input"]
    database = arg_dict["database"]
    perc_path = arg_dict["percolator_path"]
    percadapter_path = arg_dict["percolator_adapter_path"]
    outfile = arg_dict["output"]

    # Run the database search, store results to idXML file
    protein_ids, peptide_ids = sse_algorithm(searchfile, database)
    sse_res_file = "sse_results.idXML"
    IdXMLFile().store(sse_res_file, protein_ids, peptide_ids)

    # TODO: Integrate predicted intensities

    # Run PercolatorAdapter
    perc_protein_ids, perc_peptide_ids = run_percolator(sse_res_file, perc_path, percadapter_path)

    # FDR filtering
    perc_peptide_ids_filtered = fdr_filtering(database, perc_peptide_ids)

    # Store result
    IdXMLFile().store(outfile, perc_protein_ids, perc_peptide_ids_filtered)


def sse_algorithm(searchfile: str, database: str):
    # Run SimpleSearchEngineAlgorithm, store protein and peptide ids

    protein_ids = []
    peptide_ids = []

    # Enable additional score annotations
    simplesearch = SimpleSearchEngineAlgorithm()
    params = simplesearch.getDefaults()

    score_annot = [b'fragment_mz_error_median_ppm', b'precursor_mz_error_ppm']
    params.setValue(b'annotate:PSM', score_annot)

    simplesearch.setParameters(params)

    simplesearch.search(searchfile, database, protein_ids, peptide_ids)

    return protein_ids, peptide_ids


def run_percolator(sse_results: str, perc_path: str, percadapter_path: str):

    # Define the command for the PercolatorAdapter run
    percadapter_command = percadapter_path + " -in " + sse_results + " -out sse_results_percolated.idXML " + \
                          "-percolator_executable " + perc_path + " -train_best_positive -score_type q-value " + \
                          "-generic_feature_set"

    # Command with pin and weight outputs
    """
    percadapter_command = percadapter_path + " -in " + sse_results + " -out sse_results_percolated.idXML " + \
                          "-percolator_executable " + perc_path + " -out_pin SSE_results_percolated.tab " + \
                          "-weights SSE_results_percolated.weights -train_best_positive -score_type q-value " + \
                          "-generic_feature_set"
    """

    os.system(percadapter_command)

    # Load the new ids
    perc_protein_ids = []
    perc_peptide_ids = []

    IdXMLFile().load("sse_results_percolated.idXML", perc_protein_ids, perc_peptide_ids)

    return perc_protein_ids, perc_peptide_ids


def fdr_filtering(database: str, peptide_ids: list):
    # Load fasta file
    fasta_file = FASTAFile()
    fasta_entries = []
    fasta_file.load(database, fasta_entries)

    # Annotate q-value
    fdr = FalseDiscoveryRate()
    fdr.apply(peptide_ids)

    # Filter by 1% PSM FDR
    idfilter = IDFilter()
    idfilter.filterHitsByScore(peptide_ids, 0.01)
    idfilter.removeDecoyHits(peptide_ids)

    return peptide_ids



if __name__ == "__main__":
    main()