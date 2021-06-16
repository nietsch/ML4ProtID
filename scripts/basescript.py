import os
from CTDopts.CTDopts import CTDModel
from CTDsupport import *
from pyopenms import *
import pandas as pd


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

    # Register file containing predicted intensities (Prosit output in generic text format)
    model.add(
        "predicted_intensities",
        required=True,
        type="input-file",
        is_list=False,
        file_formats=["generic"],
        description="Generic text file containing predicted peak intensities (Prosit output)"
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

    # Register output file name (FDR filtered idXML file)
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
    predicted_intensities = arg_dict["predicted_intensities"]
    perc_path = arg_dict["percolator_path"]
    percadapter_path = arg_dict["percolator_adapter_path"]
    outfile = arg_dict["output"]

    # Run the database search, store results to idXML file
    protein_ids, peptide_ids = sse_algorithm(searchfile, database)
    sse_res_file = "sse_results.idXML"
    IdXMLFile().store(sse_res_file, protein_ids, peptide_ids)

    # Generate theoretical spectra for the (filtered) hits found by database search
    peptide_ids_filtered = filter_peptides(peptide_ids)
    theoretical_exp = theoretical_spectra(peptide_ids_filtered)

    # Integrate predicted intensities to theoretical spectra
    theoretical_exp_intensities = integrate_intensities(predicted_intensities, theoretical_exp)

    # Align experimental and theoretical spectra
    experimental_exp = MSExperiment()
    MzMLFile().load(searchfile, experimental_exp)
    spectrum_alignment(experimental_exp, theoretical_exp_intensities, peptide_ids_filtered)

    # TODO: Compute and add new meta value to hits

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


def filter_peptides(peptide_ids: list):
    # Filter peptides according to peptide sequence length limit of Prosit
    peptide_ids_filtered = []

    for p in peptide_ids:
        hits = []
        for h in p.getHits():
            sequence = str(h.getSequence())

            # Skip sequences exceeding the limit of 30 amino acids
            if len(h.getSequence().toUnmodifiedString()) > 30:
                continue

            hits.append(h)

        p.setHits(hits)
        if hits:
            peptide_ids_filtered.append(p)

    return peptide_ids_filtered


def theoretical_spectra(peptide_ids: list):
    # Generate theoretical spectra

    tsg = TheoreticalSpectrumGenerator()
    theoretical_exp = MSExperiment()

    for p in peptide_ids:
        for hit in p.getHits():
            spec = MSSpectrum()
            sequence = str(hit.getSequence())

            peptide = AASequence.fromString(sequence)

            p = Param()
            p.setValue("add_b_ions", "true")
            p.setValue("add_metainfo", "true")

            tsg.setParameters(p)

            # Generate ions with a maximum charge of 3
            tsg.getSpectrum(spec, peptide, 1, min(hit.getCharge(), 3))

            theoretical_exp.addSpectrum(spec)

    return theoretical_exp


def integrate_intensities(generic_out: str, theoretical_exp: MSExperiment):
    # Parse Prosit output (given in generic text format)
    # TODO: Work with dictionary instead of list for faster lookup

    df = pd.read_csv(generic_out)

    # Get all rows (i.e. ions) associated with a hit and store them as an element in a list in order to have a better
    # access to needed values in the following mapping step
    predicted_peaks = []
    tmpspec = []
    tmp = df.at[0, 'PrecursorMz']

    for i, r in df.iterrows():

        if r['PrecursorMz'] != tmp:
            predicted_peaks.append(tmpspec)
            tmpspec = []
            tmp = df.at[i, 'PrecursorMz']
            continue

        tmpspec.append(r)
        tmp = r['PrecursorMz']

    predicted_peaks.append(tmpspec)

    # Map the predicted intensities to the respective peaks by matching the fragment types and numbers
    # Store the adjusted spectra in a new experiment
    ints_added_exp = MSExperiment()
    idx = 0

    for s in theoretical_exp:

        pred_mz = []
        pred_int = []

        for ion, peak in zip(s.getStringDataArrays()[0], s):
            for r in predicted_peaks[idx]:
                if ion.decode()[0] == r['FragmentType'] \
                        and int(ion.decode()[1:].split('+', 1)[0]) == r['FragmentNumber'] \
                        and ion.decode().count('+') == r['FragmentCharge']:

                    pred_mz.append(peak.getMZ())
                    pred_int.append(r['RelativeIntensity'])

        s.set_peaks((pred_mz, pred_int))
        ints_added_exp.addSpectrum(s)

        idx += 1

    return ints_added_exp


def spectrum_alignment(experimental_exp: MSExperiment, theoretical_exp_intensities: MSExperiment, peptide_ids: list):
    # Align experimental and theoretical spectra

    # Create dictionary for assigning indices to spectra native IDs
    spectrum_index = 0
    native_id2spectrum_index = dict()

    for s in experimental_exp.getSpectra():
        native_id2spectrum_index[s.getNativeID()] = spectrum_index
        spectrum_index += 1

    # Match experimental spectra and (filtered) db search results in order to allow spectrum alignment
    for pep_idx, pep in enumerate(peptide_ids):

        ident_native_id = pep.getMetaValue("spectrum_reference")
        spectrum_index = native_id2spectrum_index[ident_native_id]

        if experimental_exp[spectrum_index].getNativeID() == ident_native_id:
            alignment = []
            spa = SpectrumAlignment()
            p = spa.getParameters()
            p.setValue(b'tolerance', 0.2)
            p.setValue(b'is_relative_tolerance', b'false')
            spa.setParameters(p)

            spec_theo = theoretical_exp_intensities[pep_idx]
            spec_exp = experimental_exp[spectrum_index]

            spa.getSpectrumAlignment(alignment, spec_theo, spec_exp)

            # Print out matched peaks
            """
            print("Matched peaks: " + str(len(alignment)))
            print("ion\ttheo. m/z\tobserved m/z")
            for theo_idx, obs_idx in alignment:
                print(spec_theo.getStringDataArrays()[0][theo_idx].decode() + "\t" +
                    str(spec_theo[theo_idx].getMZ()) + "\t" +
                    str(spec_exp[obs_idx].getMZ()))
            """

            # TODO: Continue


def run_percolator(sse_results: str, perc_path: str, percadapter_path: str):

    # Define the command for the PercolatorAdapter run
    percadapter_command = percadapter_path + " -in " + sse_results + " -out sse_results_percolated.idXML " + \
                          "-percolator_executable " + perc_path + " -train_best_positive -score_type q-value " + \
                          "-generic_feature_set"

    # Command with pin and weight outputs
    """
    percadapter_command = percadapter_path + " -in " + sse_results + " -out sse_results_percolated.idXML " + \
                          "-percolator_executable " + perc_path + " -out_pin sse_results_percolated.tab " + \
                          "-weights sse_results_percolated.weights -train_best_positive -score_type q-value " + \
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