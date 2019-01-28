SMURF software package ![GitHub Logo](smurf-logo.jpg)
========================
**SMURF: Short MUltiple Regions Framework** – A software package for microbial profiling based on combining sequencing results of any number of amplified 16S rRNA regions.

Installation 
------------
Download the SMURF package, join and decompress all the database files and folders as follows:
"cat
./Green_Genes_201305/unique_up_to_3_ambiguous_16S/Green_Genes_201305_unique_up_to_3_ambiguous_16S.fasta.gz*> ./Green_Genes_201305/unique_up_to_3_ambiguous_16S/Green_Genes_201305_unique_up_to_3_ambiguous_16S.fasta.gz" 
“gunzip ./Green_Genes_201305/unique_up_to_3_ambiguous_16S/Green_Genes_201305_unique_up_to_3_ambiguous_16S.fasta.gz”

Inputs
------
**Sequenced reads** – fastq files of the sequenced reads. An example sample that was sequenced on paired-end Illumina MiSeq is provided.  

**Ad hoc database** – a database of k-mers per region (an example of the DB for the six primer pairs described in the manuscript is provided with the software package)

Output
------
The profiling results is a csv file with reconstructed groups' information. A group is a set of full-length 16S rRNA sequences that share their sequence over the amplified regions.
Results include:

**Group frequency** – the frequency assigned to each group

**Read count** – the number of reads assigned to each group

**Number of sequences** – the number of sequences that belong to each group (i.e., the ambiguity)

**Taxonomy** – Classification based on the taxonomy file provided within the package. Alternatively, when using a custom DB without taxonomy file a directory "groups" is generated which contains one fasta file for each reconstructed group. The fasta file contains all the 16S sequences belonging to a group.

Usage
-------
1. Prepare the unique sequences database by running "DB_prepare.m" script. The parameters that define the DB are specified in "db_params_script.m" file
2. Generate the ad-hoc database of k-mers by running the "adhoc_DB_prepare.m script". The parameters that define the DB are specified in "adhoc_db_params_script.m" file
3. To profile a single sample use a script named "profile_one_sample.m". Parameters of the sample are specified in the script. The sample preparation and the algorithm parameters are specified in the script called params_script.m 
4. Alternatively you can run the main_smurf.m function with one config file "compiled_config_script.m" as described for the standalone versions.

Standalone versions
--------------------
To allow using SMURF when Matlab is not available, we prepared a Windows and a Linux precompiled versions of SMURF. 
The user should first download and install one of the following file:

**Windows** - download and install MCR from: https://www.mathworks.com/supportfiles/downloads/R2015a/deployment_files/R2015a/installers/win64/MCR_R2015a_win64_installer.exe

**Linux** - download and install MCR from: https://www.mathworks.com/supportfiles/downloads/R2014a/deployment_files/R2014a/installers/glnxa64/MCR_R2014a_glnxa64_installer.zip

The SMURF standalone version has a single input argument: the name of the configuration file, to be provided in the command line. 
In Windows the user will be prompted to select the config file in case it was not provided in the command line. 

An example of a configuration file (compiled_config_script.m) is provided with the package.
The software will write all its log file in the same directory.

For example, the Linux command line would be:
./SMURF_lin ../Configs/compiled_params_script.m

And for Windows:
SMURF_win.exe ../Configs/compiled_params_script.m  


Database parameters
----------------------
**primers_seq** - the list of (possibly degenerate) primers

**DB_kmer_len** - the length of the kmers saved in the ad hoc data base. Notice that the **kmer_len** parameter used for the sample profiling will have to be smaller than **DB_kmer_len**

**allowed_mm** - the maximal allowed mismatch between a primer and 16S sequence for the latter to be considered amplified by the primer 

Sample parameters
----------------------

**base_samples_dir**  - the samples’ directory 

**sample_name** – the directory of fastq files of the specific samples. Notice that the fastq files must be named using the following convention. If the sample_name=’Example’, then for paired end sequencing the files will be names: Example_L001_R1_001.fastq and Example_L001_R2_001.fastq

**primer_set_name** – the name of the primers set used for sequencing the sample

**kmer_len** – the length of the k-mer to be used for profiling

Data preprocessing parameters
-----------------------------
**data_type** – specify whether the reads contain quality scores. Possible values are ‘fasta’ or ‘fastq’.

**pe_flag** - 0/1 flag specifying whether sequencing was single-end or paired-end, respectively.

**qual_th** – minimal Phred quality to be used in "prc_high_qual".

**prc_high_qual** – minimal required proportion of nucleotides in a read, having a Phred score above qual_th. 

**low10_th** – maximal number of base pairs allowed to have Phred score below 10.

**max_num_Ns** – maximal number of ambiguous nucleotides allowed per read. 

**algo_pe_flag** – 0/1 flag specifying whether the reconstruction should be performed assuming single-end or 
paired-end, respectively. This parameter will usually be equal to pe_flag, although pe_flag = 1 and algo_pe_flag = 0 is allowed, while pe_flag = 0 and algo_pe_flag = 1 is not possible.

**max_err_inprimer** - maximal number of mismatches allowed between the primer sequence and the read. Used for assigning a read to a region.

**with_primer_flag** – 0/1 flag specifying whether to remove the primers after assignment of read to regions.

Algorithm parameters
--------------------
**uniS16_dir** – directory of the reference database used for profiling.

**db_filename** – name of the reference (fasta) database file (without extension). 

**filter_reads** – 0/1 flag specifying whether to apply the low abundance data preprocessing filter.

**min_read_freq** - minimal required frequency for a read per region to pass the low frequency reads filter.

**min_read_count** - minimal required count for a read per region to pass the low frequency reads filter.

**nMM_cut** – maximal number of mismatches allowed when matching reads to k-mers. 

**do_filter** - 0/1 flag specifying whether to apply the data preprocessing bacteria filter.

**pe** – probability of error per nucleotide assumed by the algorithm.

**tol** – maximal L1 change in the estimate of read proportions vector between algorithm iterations.

**numIter** – maximal number of iteration of the reconstruction algorithm.

