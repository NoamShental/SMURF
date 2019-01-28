function reads_preparation_func(fname, sample_num, primers_seq, PrepConfig, readsStatsObj)

% Extract params
algo_pe_flag = PrepConfig.algo_pe_flag;
max_err_inprimer = PrepConfig.max_err_inprimer;
with_primer_flag = PrepConfig.with_primer_flag;

% Degenerate primers
for rr = 1:size(primers_seq,1)
    for ii =1:2
        maxNonACGT = sum(primers_seq{rr,ii}~='A' & primers_seq{rr,ii}~='C' & primers_seq{rr,ii}~='G' & primers_seq{rr,ii}~='T');
        primers_seq{rr,ii} = unambiguit_one_seq(primers_seq{rr,ii},'',maxNonACGT);
    end
end

if PrepConfig.pe_flag == 0
    algo_pe_flag = 0;
end


% Quality filter
disp('Doing quality filters')
read_fastq_save_unireads(fname, PrepConfig, algo_pe_flag, readsStatsObj)

% Split the reads to regions
read_unireads_save_split_to_regions(fname,sample_num, primers_seq, max_err_inprimer, with_primer_flag, algo_pe_flag, readsStatsObj)



