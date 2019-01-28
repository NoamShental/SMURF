% ********************** GENERAL PARAMETERS ********************
base_samples_dir = 'C:\Users\Gari\Documents\WIS\Bacteria\Algo Paper\SMURF code';
sample_name = 'Example';
kmer_len = 100;


% ********************** SAMPLE PREP PARAMETERS ********************
% Set the 16S reference DB
uniS16_dir = 'C:/Users/Gari/Documents/WIS/Bacteria/Algo Paper/SMURF code/Green_Genes_201305/unique_up_to_3_ambiguous_16S';
db_filename = 'GreenGenes_201305_unique_up_to_3_ambiguous_16S';

% The primers seqs - 6 regions primers
primer_set_name = 'amp6Regions';
primers_seq{1,1} = 'TGGCGGACGGGTGAGTAA';
primers_seq{1,2} = 'CTGCTGCCTCCCGTAGGA';
primers_seq{2,1} = 'TCCTACGGGAGGCAGCAG';
primers_seq{2,2} = 'TATTACCGCGGCTGCTGG';
primers_seq{3,1} = 'CAGCAGCCGCGGTAATAC';
primers_seq{3,2} = 'CGCATTTCACCGCTACAC';
primers_seq{4,1} = 'AGGATTAGATACCCTGGT';
primers_seq{4,2} = 'GAATTAAACCACATGCTC';
primers_seq{5,1} = 'GCACAAGCGGTGGAGCAT';
primers_seq{5,2} = 'CGCTCGTTGCGGGACTTA';
primers_seq{6,1} = 'AGGAAGGTGGGGATGACG';
primers_seq{6,2} = 'CCCGGGAACGTATTCACC';

DB_kmer_len = 130;

% Generate kmers DB path and filename
allowed_mm = 2;
suffix = ['_' primer_set_name '_' num2str(allowed_mm) 'mm_RL' num2str(DB_kmer_len)];


% ********************** SAMPLE PREP PARAMETERS ********************
% Quality filter
PrepConfig.data_type = 'fastq';
PrepConfig.qual_th = 30;
PrepConfig.prc_high_qual = 0.75; 
PrepConfig.low10_th = 3;
PrepConfig.read_len = kmer_len;
PrepConfig.pe_flag = 1;
PrepConfig.max_num_Ns = 0; 

% Algo related
PrepConfig.algo_pe_flag = 1;
PrepConfig.max_err_inprimer = 2;
PrepConfig.with_primer_flag = 0;


% ********************** ALGORITHM  PARAMETERS ********************
AlgoConfig.verbose = 1;

% DB path
AlgoConfig.dbPath = [uniS16_dir suffix];
AlgoConfig.dbFileName = [db_filename suffix];
AlgoConfig.taxaFileName = 'gg_rdp_taxa_name_calls.mat';

% Reads filter params
AlgoConfig.filter_reads = 1;
AlgoConfig.min_read_freq = 1e-4;
AlgoConfig.min_read_count = 2;

% Reads to kmers alignment params
AlgoConfig.pe = 0.005;
AlgoConfig.nMM_cut = 2; 


% Reconstruction params
AlgoConfig.do_filter = 1;

AlgoConfig.read_type = char('PE'*PrepConfig.algo_pe_flag + 'SE'*(1-PrepConfig.algo_pe_flag));
AlgoConfig.readLen = kmer_len;
AlgoConfig.barcoded_regions = 1;

AlgoConfig.tol = 5e-7;
AlgoConfig.numIter = 10000;


