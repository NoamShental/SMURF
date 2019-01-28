
clear

% ************************* Sample's PARAMETERS *********************
base_samples_dir = './..';
sample_name = 'Example';
primer_set_name = 'amp6Regions';
kmer_len = 100;

% *********************** LOAD METHOD's PARAMETERS *******************
run('../Configs/params_script')


% *************************** LOAD DB and TAXONOMY FILE *********************
% Load the taxonomy
if exist([uniS16_dir '/' AlgoConfig.taxaFileName],'file')
    Taxonomy = load([uniS16_dir '/' AlgoConfig.taxaFileName]);
else
    Taxonomy = [];
end


% Load the 16S sequences
if isempty(Taxonomy) && ~exist('Sequence_uni','var')
    uniS16_file = [uniS16_dir '/' db_filename '.fasta'];
    [Header_uni, Sequence_uni] = fastaread(uniS16_file);
else
    Header_uni = {};
    Sequence_uni = {};
end

% *************************** PROFILE ONE SAMPLE *********************
SampleConfig = struct;
SampleConfig.sample_name = sample_name;
SampleConfig.sample_dir = [base_samples_dir '/' sample_name];
SampleConfig.primers_seq = primers_seq;

main_multiple_regions(PrepConfig,AlgoConfig,SampleConfig,Header_uni,Sequence_uni,Taxonomy)
