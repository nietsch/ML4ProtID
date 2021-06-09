from pyopenms import *
import pandas as pd

# Script for generating theoretical spectra including predicted intensities

# TODO: Integrate into basescript

# Load database search results
peptide_ids = []
IdXMLFile().load('../files/sse_results.idXML', [], peptide_ids)

# Load experimental spectra
experimental_exp = MSExperiment()
MzMLFile().load('/peptideatlas/ftp.peptideatlas.org/distro/mzML/JD_06232014_sample1-A.mzML', experimental_exp)


# Generate theoretical spectra for the hits found by the previous database search
tsg = TheoreticalSpectrumGenerator()
theoretical_exp = MSExperiment()

for p in peptide_ids:
    for hit in p.getHits():
        spec = MSSpectrum()
        sequence = str(hit.getSequence())

        # Remove (Carbamidomethyl) notation after C's, since each C is treated as C with carbamidomethylation in Prosit
        sequence = sequence.replace("(Carbamidomethyl)", "")

        # Ensure accordance with Prosit output (sequences with length > 30 are excluded)
        if len(sequence.replace("Oxidation", "ox")) <= 30:
            peptide = AASequence.fromString(sequence)

            p = Param()
            p.setValue("add_b_ions", "true")
            p.setValue("add_metainfo", "true")

            tsg.setParameters(p)
            tsg.getSpectrum(spec, peptide, 1, 1)

            print("Spectrum 1 of", peptide, "has", spec.size(), "peaks.")

            theoretical_exp.addSpectrum(spec)

MzMLFile().store("../files/theoretical_spectra.mzML", theoretical_exp)


# Parse Prosit output (given in generic text format)
df = pd.read_csv('../files/prositfiles/testrun.generic')

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
        newspec = MSSpectrum()
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
                    and int(ion.decode()[1:len(ion.decode())-1]) == r['FragmentNumber'] \
                    and r['FragmentCharge'] == 1:

                pred_mz.append(peak.getMZ())
                pred_int.append(r['RelativeIntensity'])

    s.set_peaks((pred_mz, pred_int))
    ints_added_exp.addSpectrum(s)

    idx += 1

MzMLFile().store("../files/theoretical_spectra_int.mzML", ints_added_exp)
