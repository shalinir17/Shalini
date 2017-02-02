#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def cutadapt_miRNA(sample, cutadapt_miRNA_flag):

    spawn_job(jobname = 'cutadapt_miRNA', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "120:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/cutadapt_drmaa.sh", args_list = [sample, p.RAW_DATA_DIR, p.ADAPTER, p.TRIMMED_DATA_PATH, p.PYTHON_VERSION])
    job_status(jobname = 'cutadapt_miRNA', resultspath = p.TRIMMED_DATA_PATH, SAMPLE = sample, outputfilename = sample + "_trimmed.fastq", FLAG_PATH = p.FLAG_PATH)
    return
if __name__ == '__main__':
    cutadapt_miRNA(sample, cutadapt_miRNA_flag)
    sys.exit(0)