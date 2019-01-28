function reconstruction_func(readsPath, sample_name, AlgoConfig, readsStatsObj)

% Clean the results directory
resDir = [readsPath '/resDir'];
if ~exist(resDir,'dir')
    mkdir(resDir)
else
    delete([resDir '/*'])
end
    

% *************************** LOAD KMERS DB ********************
if ~isfield(AlgoConfig,'const_len_flag')
    AlgoConfig.const_len_flag = strcmp(AlgoConfig.read_type,'SE');
end
bactData = load_bact_DB(AlgoConfig.dbPath,AlgoConfig.dbFileName,AlgoConfig.with_primer_flag,AlgoConfig.primers_len,AlgoConfig.readLen,AlgoConfig.const_len_flag,AlgoConfig.use_regions);



% *************************** LOAD EXPERIMENTAL READS ********************
files = dir([readsPath '/sample_' sample_name '_region*.mat']);
nR = length(AlgoConfig.use_regions);
for rr = 1:nR
    tmp = load([readsPath '/sample_' sample_name '_region_' num2str(AlgoConfig.use_regions(rr)) '_unireads.mat']);
    
    % Make sure the reads are unique
    [S,xi] = sortrows(tmp.readsuni);
    freq1 = tmp.frequni(xi);
    [Suni, ia] = unique(S,'rows', 'first');
    [~, ib] = unique(S,'rows', 'last');
    cumcount = [0;cumsum(freq1)];
    freq = cumcount(ib+1)-cumcount(ia);
    
    
    experimental_reads(rr).uniqueReads = Suni;
    experimental_reads(rr).uniqueReads_count = freq;
end


% *************************** RUN THE ALGORITHM ********************
% Build matrices Ad
dat0 = build_A_matrices(bactData, experimental_reads, AlgoConfig,readsStatsObj);
dat0_filename = [resDir '/sample_' sample_name '_dat0.mat'];
save(dat0_filename, 'dat0')

% Solve the mixture
[bact_freq_vec, bactMetaGroups, indices] = solve_iterative_noisy(dat0, AlgoConfig,readsStatsObj);
found_bacteria.frequency = bact_freq_vec';


% Calculate the number of reads assigned to each bacteria
read_count = zeros(length(indices),nR);
for rr = 1:nR
    if ~isempty(dat0.A{rr})
        Pr_r_and_j = bsxfun(@times,dat0.A{rr}(:,indices),found_bacteria.frequency);
        Pr_j_given_r = bsxfun(@rdivide,Pr_r_and_j,eps+sum(Pr_r_and_j,2));
        count_j_given_r = bsxfun(@times,Pr_j_given_r,dat0.F{rr});
        read_count(:,rr) = round(sum(count_j_given_r,1)');
    end
end
found_bacteria.assigned_reads = sum(read_count,2)';

% Save results
matlab_filename = [resDir '/sample_' sample_name '_results.mat'];
save(matlab_filename, 'AlgoConfig','bactMetaGroups','found_bacteria')
%---------------------------------------------------------------------------------


