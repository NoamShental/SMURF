% Run the DB config
run('../Configs/adhoc_db_params_script')

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


