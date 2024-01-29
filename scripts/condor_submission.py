import os
import string
import numpy as np

"""
Submission scripts handle submitting analysis and pdf making scripts to condor.
"""

def submit_analysis(isotope, fv_cut, z_cut):
    # define the inputs and the outputs
    input_fpath  = "/data/snoplus3/SNOplusData/production/multisite_directionality_0p6gL_RAT_6.18.9/ratds"
    output_fpath = "/data/snoplus3/hunt-stokes/multisite_directionality_0p6gL_refactored/run_by_run_test"

    # path to save the individual sh and submit files to
    condor_path = "/data/snoplus3/hunt-stokes/multisite_directionality_0p6gL_refactored/condor"
    
    # load up the pdf run list
    run_list = np.loadtxt("../runlists/test_runlist.txt", dtype = int)

    # load up the templates
    with open("../condor/analysis_template.sh", "r") as infile:
        rawTextSh = string.Template(infile.read())
    with open("../condor/analysis_template.submit", "r") as infile:
        rawTextSubmit = string.Template(infile.read())

    for irun in run_list:
        outTextSh = rawTextSh.substitute(ISOTOPE = isotope, RUN_NUMBER = irun, FV_CUT = fv_cut, Z_CUT = z_cut, INPUT = input_fpath, OUTPUT = output_fpath)
        
        name = f"analysis_{isotope}_{irun}"
         
        outTextSubmit = rawTextSubmit.substitute(SH_NAME = f"{condor_path}/sh/{name}.sh", LOG_NAME = name)

        # create the files
        with open(f"{condor_path}/sh/{name}.sh", "w") as outfile:
            outfile.write(outTextSh)
        os.chmod(f"{condor_path}/sh/{name}.sh", 0o0777)

        with open(f"{condor_path}/submit/{name}.submit", "w") as outfile:
            outfile.write(outTextSubmit)
        
        # submit this pdf extraction
        command = f"condor_submit -b analysis_{isotope} {condor_path}/submit/{name}.submit"
        os.system(command)

def submit_pdf_maker(isotope, fv_cut, z_cut):
    """
    PDF Files are make on a run_by_run basis. So, create individual .sh and .submit
    files for each run in the PDF MC list.
    """

    # define the inputs and the outputs
    input_fpath  = "/data/snoplus3/SNOplusData/production/multisite_directionality_0p6gL_RAT_6.18.9/ratds"
    output_fpath = "/data/snoplus3/hunt-stokes/multisite_directionality_0p6gL_refactored/run_by_run_pdf"

    # path to save the individual sh and submit files to
    condor_path = "/data/snoplus3/hunt-stokes/multisite_directionality_0p6gL_refactored/condor"
    
    # load up the pdf run list
    run_list = np.loadtxt("../runlists/pdf_runlist.txt", dtype = int)

    # load up the templates
    with open("../condor/pdf_maker_template.sh", "r") as infile:
        rawTextSh = string.Template(infile.read())
    with open("../condor/pdf_maker_template.submit", "r") as infile:
        rawTextSubmit = string.Template(infile.read())

    for irun in run_list:
        outTextSh = rawTextSh.substitute(ISOTOPE = isotope, RUN_NUMBER = irun, FV_CUT = fv_cut, Z_CUT = z_cut, INPUT = input_fpath, OUTPUT = output_fpath)
        
        name = f"pdf_maker_{isotope}_{irun}"
         
        outTextSubmit = rawTextSubmit.substitute(SH_NAME = f"{condor_path}/sh/{name}.sh", LOG_NAME = name)

        # create the files
        with open(f"{condor_path}/sh/{name}.sh", "w") as outfile:
            outfile.write(outTextSh)
        os.chmod(f"{condor_path}/sh/{name}.sh", 0o0777)

        with open(f"{condor_path}/submit/{name}.submit", "w") as outfile:
            outfile.write(outTextSubmit)
        
        # submit this pdf extraction
        command = f"condor_submit -b pdf_maker_{isotope} {condor_path}/submit/{name}.submit"
        os.system(command)

# submit_pdf_maker("Tl208", 5000.0, 1000.0)
# submit_analysis("B8_Solar_Nue", 5000.0, 1000.0)
submit_analysis("Tl208", 5000.0, 1000.0)