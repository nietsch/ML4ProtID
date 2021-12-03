# Already running Prosit server needed
# Example command:
# python basescript_MSGF_RT.py -input "searchfile.mzML" -database "TargDecoy.fasta"
# -prosit_server_ip "http://x.x.x.x:xxxx" -openmstools_path ".../" -msgfplus_path ".../MSGFPlus.jar"
# -percolator_path ".../percolator" -output "output.idXML"

import os
from CTDopts.CTDopts import CTDModel
from CTDsupport import *
from pyopenms import *
import pandas as pd
import csv
from deeplc import DeepLC

def main():

    # register command line arguments
    model = CTDModel(
        name="ML4ProtID",
        version="1.0",
        description="Script for finding peptide candidates incorporating predicted theoretical peak intensities and"
                    "retention times (Spectral angle, RT difference and absolute pred. RT as additional meta values). "
                    "Database search using MSGF+ (via MSGFPlusAdapter). Rescoring with Percolator.",
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

    # Register path to MSGF+ executable, needed by OpenMS tool MSGFPlusAdapter
    model.add(
        "msgfplus_path",
        required=True,
        type="string",
        is_list=False,
        description="Path to MSGF+ executable"
    )

    # Register path to directory including OpenMS tool executables in order to run the tools with os.system
    model.add(
        "openmstools_path",
        required=True,
        type="string",
        is_list=False,
        description="Path to OpenMS tool executables directory"
    )

    # Register number of top hits to keep for each spectrum
    model.add(
        "top_hits_per_spectrum",
        required=False,
        is_list=False,
        default=1,
        description="Number of top hits to keep for each spectrum. Default: 1."
    )

    # Register normalized collision energy considered in the peak intensity prediction with Prosit
    model.add(
        "ce",
        required=False,
        is_list=False,
        default=27,
        description="Collision energy considered in peak intensity prediction (Prosit). Default: 27"
    )

    # Register IP of already running Prosit server
    model.add(
        "prosit_server_ip",
        required=True,
        type="string",
        is_list=False,
        description="Prosit server IP in the format http://x.x.x.x:xxxx"
    )

    # Register path to Percolator executable, needed by OpenMS tool PercolatorAdapter
    model.add(
        "percolator_path",
        required=True,
        type="string",
        is_list=False,
        description="Path to Percolator executable"
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
    msgfplus_path = arg_dict["msgfplus_path"]
    openmstools_path = arg_dict["openmstools_path"]
    hits_per_spec = str(arg_dict["top_hits_per_spectrum"])
    ce = int(arg_dict["ce"])
    server_ip = arg_dict["prosit_server_ip"]
    perc_path = arg_dict["percolator_path"]
    outfile = arg_dict["output"]

    # Run the database search on experimental spectra, store results to idXML file
    protein_ids, peptide_ids = msgf_algorithm(searchfile, database, msgfplus_path, openmstools_path, hits_per_spec)

    # Generate input csv file for Prosit and start the run with already running Prosit server
    generate_csv_file(peptide_ids, ce)
    prosit_command = "curl -F \"peptides=@prosit_input.csv\" " + server_ip + "/predict/generic > " \
                     "pred_ints.generic"
    os.system(prosit_command)

    # Generate calibrants, predict absolute RT and add as new meta value for all peptides
    peptide_ids_calinput = peptide_ids[:]
    calibrants = create_calibration_data(peptide_ids_calinput)
    prot_ids, pep_ids = annotate_predictions(protein_ids, peptide_ids, predict(peptide_ids_to_dataframe(peptide_ids),
                                                                               calibrants))

    # Generate theoretical spectra for the hits found by database search
    theoretical_exp, peptide_seqs = theoretical_spectra(pep_ids)

    # Integrate predicted intensities to theoretical spectra
    theoretical_exp_intensities = integrate_intensities("pred_ints.generic", theoretical_exp, peptide_seqs)

    # Align experimental and theoretical spectra, add spectral angle, RT difference and absolute predicted RT as
    # additional meta values
    experimental_exp = MSExperiment()
    MzMLFile().load(searchfile, experimental_exp)
    peptide_ids_add_vals = spectrum_alignment(experimental_exp, theoretical_exp_intensities, prot_ids, pep_ids)
    msgf_res_add_vals_file = "msgf_results_add_vals.idXML"
    IdXMLFile().store(msgf_res_add_vals_file, prot_ids, peptide_ids_add_vals)

    # Add MSGF+ specific features and run PercolatorAdapter
    perc_protein_ids, perc_peptide_ids = run_percolator(msgf_res_add_vals_file, perc_path, openmstools_path)

    # FDR filtering
    perc_peptide_ids_filtered = fdr_filtering(perc_peptide_ids)

    # Store result
    IdXMLFile().store(outfile, perc_protein_ids, perc_peptide_ids_filtered)


def msgf_algorithm(searchfile: str, database: str, msgfplus_path: str, openmstools_path: str, hits_per_spec: str):
    # Run MSGF+ using MSGFPlusAdapter, store protein and peptide ids

    msgf_outfile = "msgf_results.idXML"
    #msgf_outfile_addedfeats = "msgf_results_addedfeats.idXML"

    msgfplusadapter_path = openmstools_path + "MSGFPlusAdapter"

    # Set and execute MSGFPlusAdapter command
    msgf_command = msgfplusadapter_path + " -in " + searchfile + " -out " + msgf_outfile + " -executable " \
                   + msgfplus_path + " -database " + database + " -min_peptide_length '7' -max_peptide_length '30'" \
                                                                " -matches_per_spec " + hits_per_spec

    os.system(msgf_command)

    protein_ids = []
    peptide_ids = []
    IdXMLFile().load(msgf_outfile, protein_ids, peptide_ids)

    # Peptide indexing
    # Load database file
    fasta_file = FASTAFile()
    fasta_filename = database
    fasta_entries = []
    fasta_file.load(fasta_filename, fasta_entries)

    # Annotate peptides with proteins
    indexer = PeptideIndexing()
    indexer_param = indexer.getParameters()
    indexer_param.setValue("tolerance", 0.05)
    indexer.setParameters(indexer_param)
    indexer.run(fasta_entries, protein_ids, peptide_ids)

    IdXMLFile().store(msgf_outfile, protein_ids, peptide_ids)

    return protein_ids, peptide_ids


def generate_csv_file(peptide_ids: list, ce: int):
    header = ['modified_sequence', 'collision_energy', 'precursor_charge']

    # Set file name
    file_name = "prosit_input.csv"

    # Set the collision energy
    collision_energy = ce

    with open(file_name, 'w') as f:
        writer = csv.writer(f)

        # Write the header
        writer.writerow(header)

        # Iterate over the hits and write the rows containing the respective values
        for p in peptide_ids:
            for h in p.getHits():
                sequence = str(h.getSequence())

                # Adjust needed notation for oxidation (other modifications are not supported)
                sequence = sequence.replace("Oxidation", "ox")

                # Remove (Carbamidomethyl) notation after cysteins, since each C is treated as C with carbamidomethylation
                sequence = sequence.replace("(Carbamidomethyl)", "")

                # Prosit error for charges > 3
                if h.getCharge() > 3:
                    continue

                row = [sequence, collision_energy, h.getCharge()]

                # Write respective row to csv file
                writer.writerow(row)


# DeepLC functions

def peptide_ids_to_dataframe(pep_ids: list) -> pd.DataFrame:
    """Parse a given list of peptide identification to a pandas DataFrame compatible with DeepLC.
    See https://github.com/compomics/DeepLC#input-files for more information.
    Args:
        pep_ids: List containing PeptideIdentification.
    Returns:
        pandas DataFrame with three columns: seq, modifications and tr.
        The returned pandas DataFrame can be used directly by DeepLC.
    """
    columns = ["seq", "modifications", "tr"]
    sequences = []
    modifications = []
    tr = []
    for pep_id in pep_ids:
        rt = pep_id.getRT()
        for hit in pep_id.getHits():
            seq = hit.getSequence()
            sequences.append(seq.toUnmodifiedString())
            tr.append(rt)
            hit_mods = []
            for pos in range(0, seq.size()):
                residue = seq.getResidue(pos)
                if residue.isModified():
                    hit_mods.append("|".join([str(pos + 1), residue.getModificationName()]))
            modifications.append("|".join(hit_mods))
    data = {
        "seq": sequences,
        "modifications": modifications,
        "tr": tr
    }
    return pd.DataFrame(data, columns=columns)


def create_calibration_data(pep_ids: list) -> pd.DataFrame:
    #Create pandas DataFrame to be used by DeepLC for calibration.
    #Args:
    #    pep_ids: List containing PeptideIdentification

    # check if target_decoy information present
    has_target_decoy = False
    for pep_id in pep_ids:
        for hit in pep_id.getHits():
            if hit.metaValueExists("target_decoy"):
                has_target_decoy = True
    if has_target_decoy:
        # annotate q-value
        fdr = FalseDiscoveryRate()
        fdr.apply(pep_ids)

        # filter by 1% PSM FDR (q < 0.01)
        idfilter = IDFilter()
        idfilter.filterHitsByScore(pep_ids, 0.01)
        idfilter.removeDecoyHits(pep_ids)

    return peptide_ids_to_dataframe(pep_ids)


def annotate_predictions(prot_ids: list, pep_ids: list, predictions: list) -> list:
    """Annotates predictions and returns proteins and peptides.
    Adds a custom meta value containing the prediction made by DeepLC.
    Args:
        prot_ids:       the protein identifications
        pep_ids:        the peptide identifications (PSMs)
        predictions:    A list predictions to be inserted into the resulting idXML file.
    """

    preds_iter = iter(predictions)
    # Insert prediction for each hit of a peptide id as meta value
    for pep_id in pep_ids:
        new_hits = []
        for hit in pep_id.getHits():
            try:
                hit.setMetaValue("prediction", str(next(preds_iter)))
            except StopIteration:
                raise SystemExit("Error: Number of predictions and peptide hits does not match.")
            new_hits.append(hit)
        pep_id.setHits(new_hits)
    return [prot_ids, pep_ids]


def predict(peptides: pd.DataFrame, calibration: pd.DataFrame = None) -> list:
    """
    Make predications based on a given peptide DataFrame and an optional
    calibration DataFrame using DeepLC.
    """
    dlc = DeepLC()
    if calibration is None:
        return dlc.make_preds(seq_df=peptides, calibrate=False)
    else:
        dlc.calibrate_preds(seq_df=calibration)
        return dlc.make_preds(seq_df=peptides)


def theoretical_spectra(peptide_ids: list):
    # Generate theoretical spectra

    theoretical_exp = MSExperiment()

    tsg = TheoreticalSpectrumGenerator()
    p = tsg.getParameters()
    p.setValue("add_first_prefix_ion", "true")
    p.setValue("add_metainfo", "true")
    tsg.setParameters(p)

    peptide_seqs = []

    for p in peptide_ids:
        for hit in p.getHits():
            spec = MSSpectrum()
            sequence = str(hit.getSequence())

            peptide = AASequence.fromString(sequence)

            # Generate ions with a maximum charge of 3
            tsg.getSpectrum(spec, peptide, 1, min(hit.getCharge()-1, 3))

            theoretical_exp.addSpectrum(spec)

            peptide_seqs.append(hit.getSequence().toUnmodifiedString())

    return theoretical_exp, peptide_seqs


def integrate_intensities(generic_out: str, theoretical_exp: MSExperiment, peptide_seqs: list):
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
    s_idx = 0
    pred_idx = 0

    for s in theoretical_exp:

        pred_mz = []
        pred_int = []

        # Retain the StringDataArrays() after set_peaks call in order to access the ion names during spectrum alignment
        s_array = s.getStringDataArrays()

        if predicted_peaks[pred_idx][0]['StrippedPeptide'] != peptide_seqs[s_idx]:
            s_idx += 1

            for ion, peak in zip(s.getStringDataArrays()[0], s):
                pred_mz.append(peak.getMZ())
                pred_int.append(peak.getIntensity())

            s.set_peaks((pred_mz, pred_int))
            s.setStringDataArrays(s_array)

            ints_added_exp.addSpectrum(s)

            continue

        for ion, peak in zip(s.getStringDataArrays()[0], s):
            for r in predicted_peaks[pred_idx]:
                if ion.decode()[0] == r['FragmentType'] \
                        and int(ion.decode()[1:].split('+', 1)[0]) == r['FragmentNumber'] \
                        and ion.decode().count('+') == r['FragmentCharge']:
                    pred_mz.append(peak.getMZ())
                    pred_int.append(r['RelativeIntensity'])

        s.set_peaks((pred_mz, pred_int))
        s.setStringDataArrays(s_array)

        ints_added_exp.addSpectrum(s)

        s_idx += 1
        pred_idx += 1

    return ints_added_exp


def spectrum_alignment(experimental_exp: MSExperiment, theoretical_exp_intensities: MSExperiment, protein_ids: list,
                       peptide_ids: list):
    # Align experimental and theoretical spectra
    # Compute and add spectral angle as new meta value

    # Base peak normalize the peak intensities in the experimental spectra
    normalizer = Normalizer()
    normalizer.filterPeakMap(experimental_exp)

    # Sqrt normalize experimental and theoretical spectra
    sqrtmower = SqrtMower()
    sqrtmower.filterPeakMap(experimental_exp)

    sqrtmower2 = SqrtMower()
    sqrtmower2.filterPeakMap(theoretical_exp_intensities)

    # Set extra features in order to add additional values for PercolatorAdapter
    search_parameters = protein_ids[0].getSearchParameters()
    search_parameters.setMetaValue(b'extra_features', b'score,spectral_angle,RT_difference,RT_predicted')
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

    theo_spec_idx = 0

    # Match experimental spectra and db search results in order to allow spectrum alignment
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

                # Set length of the vectors
                peptide_len = len(hit.getSequence().toUnmodifiedString()) - 1
                len_vector = peptide_len * min(hit.getCharge() - 1, 3) * 2

                # Initialize vectors with zeros
                v_theo = np.zeros(len_vector)
                v_exp = np.zeros(len_vector)

                arr_lst = []
                for a in spec_theo.getStringDataArrays()[0]:
                    arr_lst.append(a)

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

                    if theo_idx >= len(arr_lst):
                        continue

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

                # Get predicted absolute RT and RT difference
                pred_RTval = float(hit.getMetaValue('prediction'))
                rt_diff = abs(pred_RTval - pep.getRT())

                # Set SA, RT difference and absolute RT as additional meta values
                hit.setMetaValue('spectral_angle', sa)
                hit.setMetaValue('RT_difference', rt_diff)
                hit.setMetaValue('RT_predicted', pred_RTval)

                new_hits.append(hit)

                theo_spec_idx += 1

        pep.setHits(new_hits)

    return peptide_ids


def run_percolator(msgf_res_add_vals_file: str, perc_path: str, openmstools_path: str):
    # Add MSGF+ specific features
    psmfeatureextractor_path = openmstools_path + "PSMFeatureExtractor"
    msgf_all_features_added = "msgf_all_features_added.idXML"

    psmfeats_command = psmfeatureextractor_path + " -in " + msgf_res_add_vals_file + " -out " + \
                       msgf_all_features_added + " -extra spectral_angle RT_difference RT_predicted"
    os.system(psmfeats_command)

    # Define the command for the PercolatorAdapter run
    percadapter_path = openmstools_path + "PercolatorAdapter"

    percadapter_command = percadapter_path + " -in " + msgf_all_features_added + \
                          " -out msgf_results_percolated.idXML -percolator_executable " + perc_path + \
                          " -out_pin msgf_results_percolated_pin.tab" + \
                          " -weights msgf_results_percolated.weights -train_best_positive -score_type q-value "

    os.system(percadapter_command)

    # Load the new ids
    perc_protein_ids = []
    perc_peptide_ids = []

    IdXMLFile().load("msgf_results_percolated.idXML", perc_protein_ids, perc_peptide_ids)

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