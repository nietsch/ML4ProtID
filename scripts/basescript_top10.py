import os
from CTDopts.CTDopts import CTDModel
from CTDsupport import *
from pyopenms import *
import pandas as pd
import csv

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

    # Run the database search on experimental spectra, store results to idXML file
    protein_ids, peptide_ids = sse_algorithm(searchfile, database)
    IdXMLFile().store("sse_results_top10.idXML", protein_ids, peptide_ids)
    # Storage of search results necessary in order to generate Prosit input (csv) file

    # Generate theoretical spectra for the hits found by database search
    theoretical_exp = theoretical_spectra(peptide_ids)

    # Integrate predicted intensities to theoretical spectra
    theoretical_exp_intensities = integrate_intensities(predicted_intensities, theoretical_exp)

    # Align experimental and theoretical spectra, add spectral angle and MSE as additional meta values
    experimental_exp = MSExperiment()
    MzMLFile().load(searchfile, experimental_exp)
    peptide_ids_add_vals = spectrum_alignment(experimental_exp, theoretical_exp_intensities, protein_ids, peptide_ids)
    sse_res_add_vals_file = "sse_results_add_vals.idXML"
    IdXMLFile().store(sse_res_add_vals_file, protein_ids, peptide_ids_add_vals)

    # Run PercolatorAdapter
    perc_protein_ids, perc_peptide_ids = run_percolator(sse_res_add_vals_file, perc_path, percadapter_path)

    # FDR filtering
    perc_peptide_ids_filtered = fdr_filtering(perc_peptide_ids)

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

    # Only include peptides of up to 30 amino acids (Prosit limit)
    params.setValue(b'peptide:max_size', 30)

    # Report (up to) top 10 scoring hits per spectrum
    params.setValue(b'report:top_hits', 10)

    simplesearch.setParameters(params)

    simplesearch.search(searchfile, database, protein_ids, peptide_ids)

    return protein_ids, peptide_ids


def theoretical_spectra(peptide_ids: list):
    # Generate theoretical spectra

    theoretical_exp = MSExperiment()

    tsg = TheoreticalSpectrumGenerator()
    p = tsg.getParameters()
    p.setValue("add_first_prefix_ion", "true")
    p.setValue("add_metainfo", "true")
    tsg.setParameters(p)

    for p in peptide_ids:
        for hit in p.getHits():
            spec = MSSpectrum()
            sequence = str(hit.getSequence())

            peptide = AASequence.fromString(sequence)

            # Generate ions with a maximum charge of 3
            tsg.getSpectrum(spec, peptide, 1, min(hit.getCharge()-1, 3))

            theoretical_exp.addSpectrum(spec)

    return theoretical_exp


def integrate_intensities(generic_out: str, theoretical_exp: MSExperiment):
    # Parse Prosit output (given in generic text format)

    df = pd.read_csv(generic_out)

    # Get all rows (i.e. ions) associated with a hit and store them as an element in a list in order to have a better
    # access to needed values in the following mapping step
    predicted_peaks = []
    tmpspec = []
    tmp_prec_mz = df.at[0, 'PrecursorMz']
    tmp_seq = df.at[0, 'ModifiedPeptide']

    for i, r in df.iterrows():

        if r['PrecursorMz'] != tmp_prec_mz or r['ModifiedPeptide'] != tmp_seq:
            predicted_peaks.append(tmpspec)
            tmpspec = []
            tmp_prec_mz = df.at[i, 'PrecursorMz']
            tmp_seq = df.at[i, 'ModifiedPeptide']
            continue

        tmpspec.append(r)
        tmp_prec_mz = r['PrecursorMz']
        tmp_seq = r['ModifiedPeptide']

    predicted_peaks.append(tmpspec)

    # Map the predicted intensities to the respective peaks by matching the fragment types and numbers
    # Store the adjusted spectra in a new experiment
    ints_added_exp = MSExperiment()
    idx = 0

    for s in theoretical_exp:

        pred_mz = []
        pred_int = []

        # Retain the StringDataArrays() after set_peaks call in order to access the ion names during spectrum alignment
        s_array = s.getStringDataArrays()

        for ion, peak in zip(s.getStringDataArrays()[0], s):
            for r in predicted_peaks[idx]:
                if ion.decode()[0] == r['FragmentType'] \
                        and int(ion.decode()[1:].split('+', 1)[0]) == r['FragmentNumber'] \
                        and ion.decode().count('+') == r['FragmentCharge']:

                    pred_mz.append(peak.getMZ())
                    pred_int.append(r['RelativeIntensity'])


        s.set_peaks((pred_mz, pred_int))
        s.setStringDataArrays(s_array)

        ints_added_exp.addSpectrum(s)

        idx += 1

    return ints_added_exp


def spectrum_alignment(experimental_exp: MSExperiment, theoretical_exp_intensities: MSExperiment, protein_ids: list,
                       peptide_ids: list):
    # Align experimental and theoretical spectra
    # Compute and add spectral angle and MSE as new meta values

    # Base peak normalize the peak intensities in the experimental spectra
    normalizer = Normalizer()
    normalizer.filterPeakMap(experimental_exp)

    # Set extra features in order to add additional values for PercolatorAdapter
    search_parameters = protein_ids[0].getSearchParameters()
    search_parameters.setMetaValue(b'extra_features', b'score,spectral_angle')
    protein_ids[0].setSearchParameters(search_parameters)

    # Create dictionary for assigning indices to spectra native IDs
    spectrum_index = 0
    native_id2spectrum_index = dict()

    for s in experimental_exp.getSpectra():
        native_id2spectrum_index[s.getNativeID()] = spectrum_index
        spectrum_index += 1

    # Create spectrum alignment
    spa = SpectrumAlignment()
    p = spa.getParameters()
    p.setValue(b'tolerance', 100.0)  # Tolerance of 100 ppm
    p.setValue(b'is_relative_tolerance', b'true')
    spa.setParameters(p)

    # Match experimental spectra and (filtered) db search results in order to allow spectrum alignment
    theo_spec_idx = 0
    for pep_idx, pep in enumerate(peptide_ids):

        ident_native_id = pep.getMetaValue("spectrum_reference")
        assert ident_native_id in native_id2spectrum_index
        spectrum_index = native_id2spectrum_index[ident_native_id]

        new_hits = []

        if experimental_exp[spectrum_index].getNativeID() == ident_native_id:

            for hit in pep.getHits():
                alignment = []
                spa = SpectrumAlignment()
                p = spa.getParameters()
                p.setValue(b'tolerance', 100.0)
                p.setValue(b'is_relative_tolerance', b'true')
                spa.setParameters(p)

                spec_theo = theoretical_exp_intensities[theo_spec_idx]
                spec_exp = experimental_exp[spectrum_index]

                spa.getSpectrumAlignment(alignment, spec_theo, spec_exp)

                # Spectral angle calculation
                peptide_len = len(hit.getSequence().toUnmodifiedString()) - 1
                len_vector = peptide_len * min(hit.getCharge() - 1, 3) * 2

                # Initialize vectors with zeros
                v_theo = np.zeros(len_vector)
                v_exp = np.zeros(len_vector)

                # Print out matched peaks
                # print("Matched peaks: " + str(len(alignment)))
                # print("ion\ttheo. m/z\ttheo. int.\tobserved m/z\t observed int.")
                for theo_idx, obs_idx in alignment:
                    # print(hit.getSequence().toUnmodifiedString(), hit.getSequence())
                    # print(spec_theo.getStringDataArrays()[0][theo_idx].decode() + "\t" +
                    #      str(spec_theo[theo_idx].getMZ()) + "\t" +
                    #      str(spec_theo[theo_idx].getIntensity()) + "\t" +
                    #      str(spec_exp[obs_idx].getMZ()) + "\t" +
                    #      str(spec_exp[obs_idx].getIntensity()))

                    # Get ion and determine index shift in vector
                    ion = spec_theo.getStringDataArrays()[0][theo_idx].decode()
                    charge_shift = ((ion.count('+') - 1) * peptide_len)

                    # Set respective index in the vectors
                    int_idx = 0
                    if ion[0] == 'b':
                        int_idx = (int(ion[1:].split('+', 1)[0]) - 1) + charge_shift
                    # Store y-ion intensities in the 2nd half of the array
                    else:
                        int_idx = (int(ion[1:].split('+', 1)[0]) - 1) + charge_shift + int(len_vector / 2)

                    # Write intensity values to the vectors
                    v_theo[int_idx] = spec_theo[theo_idx].getIntensity()
                    v_exp[int_idx] = spec_exp[obs_idx].getIntensity()

                # Normalize vectors with L2 norm
                # Avoid division by zero problems by checking if arrays are empty
                if v_theo.any():
                    v_theo_norm = v_theo / np.linalg.norm(v_theo)
                else:
                    v_theo_norm = v_theo

                if v_exp.any():
                    v_exp_norm = v_exp / np.linalg.norm(v_exp)
                else:
                    v_exp_norm = v_exp

                # Compute dot product
                v_dot = np.dot(v_theo_norm, v_exp_norm)

                # Compute inverse cosine
                v_inv_cos = np.arccos(v_dot)

                # Compute spectral angle
                sa = 1 - 2 * (v_inv_cos / np.pi)

                # Compute mean squared error (Omitted due to lack of improvement)
                #squared_diff = [(theo_int - exp_int) ** 2 for theo_int, exp_int in zip(v_theo, v_exp)]
                #mse = 1.0 # If there are no matching peaks
                #if len(alignment) > 0:
                #    mse = sum(squared_diff) / len(alignment)

                hit.setMetaValue('spectral_angle', sa)
                # hit.setMetaValue('mse', mse)

                new_hits.append(hit)

                theo_spec_idx += 1

        pep.setHits(new_hits)

    return peptide_ids


def run_percolator(sse_results: str, perc_path: str, percadapter_path: str):

    # Define the command for the PercolatorAdapter run
    percadapter_command = percadapter_path + " -in " + sse_results + " -out sse_results_percolated.idXML " + \
                          "-percolator_executable " + perc_path + " -out_pin sse_results_percolated.tab " + \
                          "-weights sse_results_percolated.weights -train_best_positive -score_type q-value "

    os.system(percadapter_command)

    # Load the new ids
    perc_protein_ids = []
    perc_peptide_ids = []

    IdXMLFile().load("sse_results_percolated.idXML", perc_protein_ids, perc_peptide_ids)

    return perc_protein_ids, perc_peptide_ids


def fdr_filtering(peptide_ids: list):
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