from pyopenms import MzMLFile, MSExperiment

def process_spectrum(spectrum, msp_file):
    precursors = spectrum.getPrecursors()

    if precursors:
        mw = precursors[0].getMZ()
        RT = spectrum.getRT()

        msp_file.write(f"Name: Spectrum\nMW: {mw}\nRT: {RT}\n")
        msp_file.write(f"Num peaks: {spectrum.size()}\n")

        for peak in spectrum:
            m_z = peak.getMZ()
            intensity = peak.getIntensity()

            try:
                annotation = peak.getMetaValue("annotation").decode("utf-8")
            except Exception as e:
                annotation = "No Annotation"

            msp_file.write(f"{m_z}\t{intensity}\t\"{annotation}\"\n")

        msp_file.write("\n")

def main():
    mzml_file_path = 'Z:/zelhamraoui/umpire-zahra/2022MQ050_ZAEL_008_01_PoolProteome_1ug_DIA_OT_CE30_Q1.mzML'
    msp_file_path = 'Z:/zelhamraoui/umpire-zahra/umpire_30CE_data_Q1.msp'

    exp = MSExperiment()
    MzMLFile().load(mzml_file_path, exp)

    with open(msp_file_path, 'w') as msp_file:
        for spectrum in exp:
            process_spectrum(spectrum, msp_file)

if __name__ == "__main__":
    main()
