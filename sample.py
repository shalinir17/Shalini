#!/usr/bin/env python
from ruffus import *
import sys 
import os
import time
import datetime 
import drmaa
###############################################
def check_file_exists(input_file, output_file):
    '''Checks if file exists'''
    if not os.path.exists(output_file):
        return True, "Missing file %s for %s" % (output_file, input_file)
    else:
        return False, "File %s exists for %s" % (output_file, input_file) 
   
def spawn_job(jobname, SAMPLE, LOG_PATH, RESULTS_EMAIL, SCHEDULER, walltime, queue, nodes, ppn, memory, script, args_list):
    '''Spawns a job on a cluster (HPC or AWS Star Cluster) using DRMAA'''
    import drmaa
    s = drmaa.Session()
    s.initialize()
    print 'Creating job template for ' + jobname
    jt = s.createJobTemplate()
    print 'Job template created'
    jt.jobName = jobname + "_" + str(SAMPLE)
    print "job name is " + jt.jobName
    jt.outputPath = LOG_PATH
    #print "Error path is" + jt.outputPath
    jt.errorPath = LOG_PATH
    jt.email = RESULTS_EMAIL
    #print "email is" + str(jt.email)
    if SCHEDULER == "PBS":
        jt.hardWallclockTimeLimit = walltime
        jt.softWallClockTimeLimit = walltime
        jt.hardRunDurationLimit = walltime
        jt.softRunDurationLimit = walltime
        jt.nativeSpecification = "-q " + queue + " -l " + "nodes=" + str(nodes) + ":ppn=" + str(ppn) + " -l mem=" + memory  
        print "native specification is" + jt.nativeSpecification
    elif SCHEDULER == "LSF":
        print "LSF currently unsupported. Please contact us to request this feature."
    elif SCHEDULER == "Slurm":
        print "Slurm currently unsupported. Please contact us to request this feature."
    elif SCHEDULER == "SGE":
        jt.nativeSpecification = "-V -cwd -pe local " + str(ppn) + " -l h_vmem=" + re.sub("gb","G",memory)
        print "native specification is " + jt.nativeSpecification
    else:
        print "Scheduler unsupported. Please make sure you have a parameter in your parameter file SCHEDULER with the options PBS, SGE, LSF, StarCluster or Slurm"

    jt.remoteCommand = os.getcwd() + script 
    #print "remote command is" + jt.remoteCommand   #THIS PRINTS THEN HANGS
    jt.args = args_list 
    print "ArgsList is " + str(jt.args)
    jt.joinFiles = True
    jobid = s.runJob(jt)
    now = datetime.datetime.now()
    date = now.strftime("%Y-%m-%d %H:%M")
    print "Date/Time: " + date 
    print "Job has been submitted with id" + jobid + " at Date/Time: " + date
    retval = s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
    now = datetime.datetime.now()
    date = now.strftime("%Y-%m-%d %H:%M")
    print "Job: " + str(retval.jobId) + ' finished with status: ' + str(retval.hasExited) + ' and exit status: ' + str(retval.exitStatus) + " at Date/Time: " + date
    print "Date/Time: " + date 
    print 'Cleaning up'
    s.deleteJobTemplate(jt)
    s.exit()
    return 

def job_status(jobname, resultspath, SAMPLE, outputfilename, FLAG_PATH):
    '''Checks to see if a job has successfully finished by checking if the specified output file exists.'''
    stage = jobname
    try:
        size = os.path.getsize(resultspath + "/" + outputfilename) 
    except OSError:
        print("Looking for file: " + resultspath + "/" + outputfilename)
        print("%s failed to produce any output files" % stage) 
    if 'size' in locals():
        if size == 0:
            print "Job Failed!"
            print('%s produced an empty output file' % stage)
        else:
            print("%s Finished and Successfully produced an output file of size %s" % (stage,size))
            flag_file = "%s/%s_%s_completed.flag" % (FLAG_PATH, stage, SAMPLE)
            open(flag_file, 'w').close()
    else:
        print("%s failed to produce any output files" % stage)   
    return


def decompress(file_directory, compression_type):
    '''Unzips (bzip or gzip) files'''
    if compression_type == "bzip":
        os.chdir(file_directory)
        os.system("bunzip2 *.bz2")
        print "Bzip files successfully uncompressed"
    elif compression_type == "gzip":
        os.chdir(file_directory)
        os.system("gunzip *.gz")
        print "Gzip files successfully uncompressed"
    else:
        print "Working with uncompressed files"
    return

def compress(file_directory, compression_type):
    '''Compresses (bzip or gzip) files'''
    if compression_type == "bzip":
        os.chdir(file_directory)
        os.system("bzip2 *.fastq")
        print "Bzip files successfully compressed"
    elif compression_type == "gzip":
        os.chdir(file_directory)
        os.system("gzip *.fastq")
        print "Gzip files successfully compressed"
    else:
        print "Working with uncompressed files" 
    return


class Bunch(object):
    '''Bunches parameters into a dictionary'''
    def __init__(self, adict):
        self.__dict__.update(adict)
        
        
def check_create_dir(directory):
    '''Creates a directory if it does not exist'''
    if not os.path.exists(directory):
        os.makedirs(directory)
    return


def parse_xml(file):
    '''Parses XML file to extract sample names from TCGA XML manifest'''
    tree = ET.parse(file)
    root = tree.getroot()
    for child in root:
        tree = ET.ElementTree(child)
        for id in tree.findall('analysis_id'):
            print id.text    
            name = "/gpfs/home/kfisch/" + id.text + ".xml"
    tree.write(name)     
    return
  

def get_TCGA_ID(file):
    '''Gets TCGA ID from TCGA XML manifest and creates a sample list'''
    tree = ET.parse(file)
    root = tree.getroot()
    sample_list = []
    for id in root.findall('Result'):
        analysis_id = id.find('analysis_id').text
        analysis_id_string = "TCGA_" + str(analysis_id)
        sample_list.append(analysis_id_string)
    print sample_list    
    return sample_list   


def make_params(step, sample_list, flag_path):
    '''Creates parameter input lists for each step in pipeline'''
    vars()['inputList_' + step] = []
    for sample in sample_list:
        vars()['inputList_' + step].append([sample, "%s/%s_%s_completed.flag" % (flag_path, step, sample)])
        #print vars()['inputList_' + step] 
    return vars()['inputList_' + step] 


def get_samples_from_txt_file(file):
    '''Creates sample list from text file of samples'''
    global SAMPLE_LIST
    sample_file = open(file, 'r')
    reader = csv.reader(sample_file)
    sample_list = [row for row in reader]
    return SAMPLE_LIST

def sumatra_start(repository, sumatra_db_path, results_path, working_dir, hg_username, sumatra_run_name, parameters):
    '''Clones the Omics Pipe repository from Bitbucket, creates a Sumatra project, and creates a Sumatra record for the current run'''
    print "sumatra_db_path is " + sumatra_db_path
    print type(sumatra_db_path)
    check_create_dir(sumatra_db_path)
    os.chdir(sumatra_db_path)
    repo1 = hgapi.Repo(repository)
    repo_path = sumatra_db_path +"/omics_pipe"
    repo= {"url":repo_path, 
           "type":"sumatra.versioncontrol._mercurial.MercurialRepository",
           "upstream":repository}
    executable= {"path":"",
                 "version": "",
                 "type":"sumatra.programs.PythonExecutable",
                 "options":"",
                 "name": "Python"}
    sumatra_launch_mode = {"working_directory": working_dir, "type": "sumatra.launch.SerialLaunchMode"}
    data_store1 = {"root":results_path, "type": "sumatra.datastore.filesystem.FileSystemDataStore"}
    database_path = sumatra_db_path + "/records/recordstore.db"
    record_store1 = {"db_file": database_path, "type": "sumatra.recordstore.django_store.DjangoRecordStore"}
    input_datastore1 = {"root": results_path, "type": "sumatra.datastore.filesystem.FileSystemDataStore"}
    while True:
        try:
            repo1.hg_clone(url = repository, path=repo_path)
            with open(repo_path + "/.hg/hgrc", "a") as myfile:
                myfile.write("[ui]\nusername= " + hg_username)         
            print "Omics pipe repository cloned to : " + repo_path
            break
        except hgapi.hgapi.HgException:
            print "Omics pipe repository already exists."
            break
    while True:
        try:
            Project(sumatra_run_name, default_repository=repo, default_executable=executable, 
                    default_launch_mode = sumatra_launch_mode, on_changed='store-diff',
                    data_store=data_store1, record_store=record_store1, input_datastore=input_datastore1)            
            print "Sumatra project created: " + sumatra_run_name + " in directory: " + sumatra_db_path
            break
        except Exception:
            print "Sumatra project already exists, loading project: " + sumatra_run_name
            break
    project = load_project(path=sumatra_db_path)
    print project
    sumatra_params = build_parameters(parameters)
    print sumatra_params
    os.chdir(repo_path)
    repo_main = "omics_pipe/main.py"
    record = project.new_record(parameters=sumatra_params, main_file=repo_main)
    print record
    return record,project

def sumatra_end(start_time, record, project):
    '''Saves the Sumatra project to the database'''
    #file1 = open("/gpfs/home/kfisch/test/test.txt", "w")
    #file1.write("test")
    #file1.close()
    record.duration = time.time() - start_time
    print record.duration
    record.output_data = record.datastore.find_new_data(record.timestamp)
    print record.output_data
    project.add_record(record)
    print project
    project.save()
    return
    ##############################################################################
def fastqc_miRNA(sample, fastqc_miRNA_flag):

    print "sample name is: ", sample 
    if p.ENDS == "PE":
        SAMPLE1 = sample + "_1"
        SAMPLE2 = sample + "_2"        
        spawn_job(jobname = 'fastqc_miRNA', SAMPLE = SAMPLE1, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "12:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/fastqc_drmaa.sh", args_list = [SAMPLE1, p.RAW_DATA_DIR,p.QC_PATH, p.FASTQC_VERSION])
        spawn_job(jobname = 'fastqc_miRNA', SAMPLE = SAMPLE2, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "12:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/fastqc_drmaa.sh", args_list = [SAMPLE2, p.RAW_DATA_DIR,p.QC_PATH, p.FASTQC_VERSION])
        job_status(jobname = 'fastqc_miRNA', resultspath = p.QC_PATH, SAMPLE = SAMPLE1,  outputfilename = SAMPLE1 + "_fastq/" + "fastqc_data.txt", FLAG_PATH = p.FLAG_PATH)
        job_status(jobname = 'fastqc_miRNA', resultspath = p.QC_PATH, SAMPLE = SAMPLE2,  outputfilename = SAMPLE2 + "_fastq/" + "fastqc_data.txt", FLAG_PATH = p.FLAG_PATH)
    else:
        spawn_job(jobname = 'fastqc_miRNA', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "12:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/fastqc_drmaa.sh", args_list = [sample, p.TRIMMED_DATA_PATH,p.QC_PATH, p.FASTQC_VERSION])
        job_status(jobname = 'fastqc_miRNA', resultspath = p.QC_PATH, SAMPLE = sample, outputfilename = sample + "_fastq/" + "fastqc_data.txt", FLAG_PATH = p.FLAG_PATH)
    
    return

if __name__ == '__main__':
    fastqc_miRNA(sample, fastqc_miRNA_flag)
    sys.exit(0)
#################################################################
default_parameters = dict(                                                                                   #MAKE SURE TO KEEP dict() and , after each variable
#PATHS
WORKING_DIR = "/home/shalini/miRNas/reference/scripts"   ,   ###PATH TO ../omics-pipeline/omics_pipe/scripts
RAW_DATA_DIR = "/gpfs/group/su/kfisch/OA/data/rawdata/NormalvsOA20130808"           ,   
RESULTS_PATH = "/gpfs/group/su/kfisch/OA/results"                                   ,
QC_PATH = "/gpfs/group/su/kfisch/OA/data/QC"                                        ,
LOG_PATH = "/gpfs/group/su/kfisch/OA/logs"                                          ,
FLAG_PATH = "/gpfs/group/su/kfisch/OA/logs"                                         ,
TOPHAT_RESULTS = "/gpfs/group/su/kfisch/OA/results/mRNA/alignments"                 ,
STAR_RESULTS = "/gpfs/group/su/kfisch/OA/results/mRNA/STAR_alignments"              ,
HTSEQ_RESULTS = "/gpfs/group/su/kfisch/OA/results/mRNA/counts"                      ,
HTSEQ_GENCODE_RESULTS = "",
CUFFLINKS_RESULTS = "/gpfs/group/su/kfisch/OA/results/mRNA/assemblies"              ,
CUFFMERGE_RESULTS = "/gpfs/group/su/kfisch/Bodymap/results/assemblies/cuffmerge"    ,
CUFFDIFF_RESULTS = "/gpfs/group/su/kfisch/Bodymap/results/cuffdiff"                 ,
TEMP_DIR = "/scratch/kfisch"                                                        ,
VARIANT_RESULTS = "/gpfs/group/su/kfisch/OA/results/variants"                       ,
FUSIONCATCHER_RESULTS = "/gpfs/group/su/kfisch/OA/results/fusions"                  ,
TOPHAT_FUSION_RESULTS = "/gpfs/group/su/kfisch/OA/results/fusions"                  ,


RESULTS_EMAIL = "kfisch@scripps.edu"                                                ,
QUEUE = "workq"                                                                   ,
SCHEDULER = "PBS",

#SAMPLE INFO
ENDS = "SE"                                                                         ,     #PE=paired-ends, SE=single ends  If paired end, sample file names must be samplename_1.fastq and samplename_2.fastq
COMPRESSION = "None",
SAMPLE_LIST =[
    'Normal_Cart_2_2',
    'Normal_Cart_3_3',
#   'Normal_Cart_4_4',
#   'Normal_Cart_5_5',
#    'Normal_Cart_6_6',
#   'OA_Cart_1_7',
#    'OA_Cart_2_8',
#    'OA_Cart_3_9',
#    'OA_Cart_4_10',
#    'OA_Cart_5_5',
#    'Normal_Cart_10_8',
#     'Normal_Cart_7_3',
#     'Normal_Cart_9_7',
#     'OA_Cart_10_9',
#     'OA_Cart_6_1',
#     'OA_Cart_7_2',
#    'OA_Cart_8_5',
#     'OA_Cart_9_6'
]                                                                 ,

#GENOME INFO
GENOME = "/gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"     ,  
REF_GENES = "/gpfs/group/databases/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"           ,
CHROM = "/gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes"                     ,
BOWTIE_INDEX = "/gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"      ,
STAR_INDEX = "/gpfs/group/databases/Homo_sapiens/UCSC/hg19/star_genome"                         ,
REF_GENES_GENCODE = "",    
BWA_INDEX = "",
#SOFTWARE VERSIONS
FASTQC_VERSION = "0.10.1"                                                           ,
TOPHAT_VERSION = "2.0.9"                                                            ,
STAR_VERSION = "2.3.0"                                                              ,
CUFFLINKS_VERSION = "2.1.1"                                                         ,            
R_VERSION = "3.0.1"                                                                 ,
SAMTOOLS_VERSION = "0.1.19"                                                         ,    
ANNOVAR_VERSION = "07292013"                                                        ,
VCFTOOLS_VERSION = "0.1.10"                                                         ,
VARSCAN_VERSION = "2.3.6"                                                           ,
FUSIONCATCHER_VERSION = "0.98"                                                      ,   
TRIMGALORE_VERSION = "0.3.3",

#SOFTWARE OPTIONS
TOPHAT_OPTIONS = "-a 5 --microexon-search --library-type fr-secondstrand --no-discordant"                                               ,    #-p8 -G options hardcoded in script. These options go after.
STAR_OPTIONS = "--readFilesCommand cat --runThreadN 8 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical"       ,  
HTSEQ_OPTIONS = "-m intersection-nonempty -s no -t exon"                                                                                ,
CUFFLINKS_OPTIONS = "-u -N"                                                                                                             ,   #'-p8 -g REF_GENES -b GENOME' options hardcoded in script

STEP = "last_function"                                                                                                                  , #enclose name in brackets (eg [last_function]) if you only want to rerun that step in the pipeline. If you want to run between certain steps in the pipeline, write steps enclosed in brackets separated by a comma (eg [STAR], [htseqcount])
STEPS = [
        'fastqc',
        'tophat', 
        'star', 
        'htseq',
        'cufflinks',
        'fusion_catcher',
        'call_variants',
        'annotate_variants',
        'last_function'
        ]                                                                                                                               ,    



PIPE_VERBOSE = 5                                                                                                                        ,
PIPE_MULTIPROCESS = 100                                                                                                                 ,   #sample number * number parallel tasks in pipeline
PIPE_REBUILD = "True"                                                                                                                   ,   #Add gnu_make_maximal_rebuild_mode = False if you want to run only one task

#Variant calling options
ANNOVARDB = "/home/shalini/miRNas/reference/humandb"                                                                                                                                     ,
VARSCAN_PATH = "/opt/applications/varscan/2.3.6/VarScan.jar"                                                                                                                            ,
VARSCAN_OPTIONS = "--min-var-freq 0.5 --min-avg-qual 30 --p-value 0.995 --output-vcf 1"                                                                                                 ,
ANNOVAR_OPTIONS = "-buildver hg19 -protocol nonsyn_splicing,1000g2012apr_all,esp6500_ea,esp6500_aa,snp135NonFlagged,cg46,ljb_sift,ljb_pp2,dominant -operation g,f,f,f,f,f,f,f,m"        ,
ANNOVAR_OPTIONS2 = "-genetype knowngene --remove"                                                                                                                                       ,
SAMTOOLS_OPTIONS = "-C 500",
#Fusion Catcher options
FUSIONCATCHERBUILD_DIR = "/gpfs/group/databases/ensembl_v72"                                                                                                                            ,
FUSIONCATCHER_OPTIONS = ""                                                                                                                                                              ,

REPORT_SCRIPT = "/gpfs/group/sanford/src/RNA/knitMeR.R",
REPORT_RESULTS = "/gpfs/group/sanford/patient/SSKT/test_patient/RNA/RNA_seq",
R_MARKUP_FILE = "/gpfs/group/sanford/src/RNA/sanfordRNASeqReport.Rmd",
BAM_FILE_NAME = "accepted_hits.bam",
R_SOURCE_PATH = "/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting",

#Trim galore options
ADAPTER = "TGGAATTCTCGGGTGCCAAGG",
TRIM_LENGTH_MIN = 10,
TRIMMED_DATA_PATH = "/gpfs/home/kfisch",

miRNA_BOWTIE_INDEX = "/gpfs/group/su/kfisch/references/bowtieindexes/miRNA",
miRNA_GTF=  "/gpfs/group/su/kfisch/references/miRBase/hsa.gff3",
CUFFMERGE_OPTIONS=  "-p 8",
CUFFMERGETOCOMPARE_OPTIONS = "-CG",
CUFFDIFF_OPTIONS=  "-p 8 -FDR 0.01 -L Normal,OA -N --compatible-hits-norm",
CUFFDIFF_INPUT_LIST_COND1= "/gpfs/group/su/kfisch/OA/results/miRNA/bamfiles/sample1.bam,/gpfs/group/su/kfisch/OA/results/miRNA/bamfiles/sample2.bam" ,
CUFFDIFF_INPUT_LIST_COND2 = "/gpfs/group/su/kfisch/OA/results/miRNA/bamfiles/sample3.bam",
#TCGA download
TCGA_XML_FILE = "/gpfs/group/sanford/TCGA/tcga_test.xml",
TCGA_KEY = "/gpfs/group/sanford/TCGA/cghub.key",
TCGA_OUTPUT_PATH = "/gpfs/group/sanford/TCGA",
SSH_USER_NAME = "kfisch@hpcdata.scripps.edu",

#BWA and SNPIR
SNPIR_RESULTS = "",
BWA_VERSION ="0.7.4",
PICARD_VERSION ="1.92",
GATK_VERSION ="3.1-1",
BEDTOOLS_VERSION ="2.17.0",
UCSC_TOOLS_VERSION ="273",
REPEAT_MASKER = "",
SNPIR_ANNOTATION ="",
RNA_EDIT  ="",
DBSNP  ="",
MILLS  ="",
G1000  ="",
BWA_RESULTS ="",
GATK_READ_GROUP_INFO ="",
SNPIR_VERSION = "1.0",
PARAMS_FILE = "",
VCF_FILE= "/gpfs/group/sanford/patient/SSKT/SSKT_2/RNA/RNA_seq/results/SNPIR/final_variants.vcf",
INTOGEN_OPTIONS= "--single-tumor ",
INTOGEN_RESULTS= "/gpfs/group/sanford/patient/SSKT/SSKT_2/RNA/RNA_seq/results/intogen",
INTOGEN_VERSION= '2.3.0',
INTOGEN_CONFIG= "/gpfs/group/sanford/patient/SSKT/SSKT_2/RNA/RNA_seq/results/intogen.conf",
USERNAME= "kfisch",
SNPIR_CONFIG = "/gpfs/group/su/meissto/cancer_report/src/RNA",
SNPIR_DIR = "/opt/applications/snpir/1.0/bin",
RSEQC_VERSION = "2.3.7",
RSEQC_REF = "/gpfs/group/su/meissto/data/knowngene.bed",
SNPEFF_VERSION = "3.5",
dbNSFP = "/gpfs/group/su/meissto/data/dbNSFP/dbNSFP2.0.txt",
SNP_FILTER_OUT_REF = "/gpfs/group/sanford/src/RNA/vcf/common_no_known_medical_impact_00-latest.vcf",
TISSUE= "EPI",
#ChIP-seq parameters
BOWTIE_VERSION = '1.0.0',
BOWTIE_OPTIONS = '',
HOMER_VERSION = '4.5',
HOMER_TRIM_OPTIONS = '-3 GATCGGAAGAGCACACGTCT -mis 1 -minMatchLength 6 -min 45',
MACS_VERSION = '1.4.2',
STEPS_PAIRS = ['macs'], 
PAIR_LIST = 'test1_chip-test1_input test2_chip-test2_input',
CHROM_SIZES = '/gpfs/group/su/kfisch/references/hg19.chrom.sizes',
MACS_RESULTS = '',
HOMER_PEAKS_OPTIONS = '-style factor -o auto' ,
HOMER_ANNOTATE_OPTIONS = '',
HOMER_MOTIFS_OPTIONS = '-size 200 -mask',
HOMER_RESULTS = "/gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/homer",

#OMICS_PIPE
DRMAA_PATH = "/opt/applications/pbs-drmaa/current/gnu/lib/libdrmaa.so",
DPS_VERSION= '1.3.1111' ,
PYTHON_VERSION="2.6.5",

#WGS parameters
BWA_OPTIONS = "-t 8 -M",
RESOURCES = "-l nodes=1:ppn=8 -l mem=31gb",

#SNIPR
TABIX_VERSION = '0.2.6',
ENCODING = "phred64",
TUMOR_TYPE= 'BRCA',
GENELIST= "/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting/ref/brca.txt",
COSMIC= "/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting/ref/cosmic.tsv",
CLINVAR="/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting/ref/clinvar.txt",
PHARMGKB_rsID= "/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting/ref/pharmgkbRSID.csv",
PHARMGKB_Allele= "/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting/ref/pharmgkbAllele.tsv",
DRUGBANK= "/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting/ref/drugbank.tsv",
CADD= "/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting/ref/cadd.tsv.gz"

)
#################################################################
p = Bunch(default_parameters)

os.chdir(p.WORKING_DIR)
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d %H:%M")    

print p

for step in p.STEPS:
    vars()['inputList_' + step] = []
    for sample in p.SAMPLE_LIST:
        vars()['inputList_' + step].append([sample, "%s/%s_%s_completed.flag" % (p.FLAG_PATH, step, sample)])
    print vars()['inputList_' + step]

@parallel(inputList_fastqc_miRNA)
@check_if_uptodate(check_file_exists)
@follows(run_fastq_length_filter_miRNA)
def run_fastqc_miRNA(sample, fastqc_miRNA_flag):
    fastqc_miRNA(sample, fastqc_miRNA_flag)
    return

def last_function(sample, last_function_flag):
    print "PIPELINE HAS FINISHED SUCCESSFULLY!!! YAY!"
    pipeline_graph_output = p.FLAG_PATH + "/pipeline_" + sample + "_" + str(date) + ".pdf"
    pipeline_printout_graph (pipeline_graph_output,'pdf', step, no_key_legend=False)
    stage = "last_function"
    flag_file = "%s/%s_%s_completed.flag" % (p.FLAG_PATH, stage, sample)
    open(flag_file, 'w').close()
    return   


if __name__ == '__main__':

    pipeline_run(p.STEP, multiprocess = p.PIPE_MULTIPROCESS, verbose = p.PIPE_VERBOSE, gnu_make_maximal_rebuild_mode = p.PIPE_REBUILD)