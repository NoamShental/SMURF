clear

% Set the primers set name (should be one of the list specified in )
primer_set_name = 'V4';

% Load the primers and pathways
run('../Configs/adhoc_db_params_script')

% Unambiguit pimers
for rr = 1:size(primers_seq,1)
    for ii =1:2
        maxNonACGT = sum(primers_seq{rr,ii}~='A' & primers_seq{rr,ii}~='C' & primers_seq{rr,ii}~='G' & primers_seq{rr,ii}~='T');
        primers_seq{rr,ii} = unambiguit_one_seq(primers_seq{rr,ii},'',maxNonACGT);
    end
end


% Load the 16S sequences
uniS16_file = [uniS16_dir '/' db_filename '.fasta'];
if ~exist('Sequence_uni','var') || ~exist('Header_uni','var') || ~exist('Dictionary','var') 
    tic
    [Header_uni,Sequence_uni] = fastaread(uniS16_file);
    toc
end
Dictionary = [];

% Create the DB dir
dbDirName = [uniS16_dir suffix];
if exist(dbDirName,'dir')
    rmdir(dbDirName, 's')
end
mkdir(dbDirName)


nB = length(Header_uni);
nR = size(primers_seq,1);

for rr = 1:nR
    
    % Generate filename
    savefilename = [dbDirName '/' db_filename suffix '_region' num2str(rr) '.mat'];
    
    % Run the insilico PCR
    save_Amat_one_region(savefilename, primers_seq(rr,:), Sequence_uni, Header_uni, Dictionary, DB_kmer_len, allowed_mm)
end


