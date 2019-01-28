function [dat0] = build_A_matrices(bactData, reads, AlgoConfig, readsStatsObj)


if ~AlgoConfig.barcoded_regions
    dat0 = [];
    return
end

pe = AlgoConfig.pe;
nMM_cut = AlgoConfig.nMM_cut;


nR = length(bactData.kmers);
nB = size(bactData.indInSeqs,1);


for rr = 1:nR
    disp(['Region ' num2str(rr) ' out of ' num2str(nR)])
    
    % Filter low abundance reads
    if AlgoConfig.filter_reads == 1
        tot_number = length(reads(rr).uniqueReads_count);
        tot_count = sum(reads(rr).uniqueReads_count);
        filter_thresh = max(AlgoConfig.min_read_freq*tot_count,AlgoConfig.min_read_count);
        remove_ind = reads(rr).uniqueReads_count < filter_thresh;
        reads(rr).uniqueReads(remove_ind,:) = [];
        reads(rr).uniqueReads_count(remove_ind) = [];
        
        addStats(readsStatsObj,'After low abundance filter', length(reads(rr).uniqueReads_count), sum(reads(rr).uniqueReads_count), rr);
        disp(['Keep high freq: ' num2str(round(100*length(reads(rr).uniqueReads_count)/tot_number)) '% of reads'])
        disp(['Keep high freq: ' num2str(round(100*sum(reads(rr).uniqueReads_count)/tot_count)) '% of counts'])
    end
    
    
    % Build M matrix
    if AlgoConfig.verbose
        disp('Building matrix M')
    end
    bamp_in_reg = find(bactData.indInSeqs(:,rr)>0);
    
    if strcmp(AlgoConfig.read_type,'PE')
        Kd = size(bactData.kmers{rr},1);
        dat0.M{rr} = sparse(bactData.indInSeqs(bamp_in_reg,rr),bamp_in_reg,ones(1,length(bamp_in_reg)),Kd,nB);
        dat0.kmers{rr} = bactData.kmers{rr};
    elseif strcmp(AlgoConfig.read_type,'SE')
        % Fwd
        fwd_readLen = AlgoConfig.readLen - (1-AlgoConfig.const_len_flag)*AlgoConfig.primers_len(rr,1);
        [values_fwd, ~, indInUni_fwd] = unique(bactData.kmers{rr}(:,1:fwd_readLen),'rows');
        Kd_fwd = length(values_fwd);
        indInValue_fwd = zeros(nB,1);
        indInValue_fwd(bamp_in_reg) = indInUni_fwd(bactData.indInSeqs(bamp_in_reg,rr));
        Ad_fwd = sparse(indInValue_fwd(bamp_in_reg),bamp_in_reg,ones(1,length(bamp_in_reg)),Kd_fwd,nB);
        
        % Rvs
        [values_rvs, ~, indInUni_rvs] = unique(bactData.kmers{rr}(:,fwd_readLen+1:end),'rows');
        Kd = length(values_rvs);
        indInValue_rvs = zeros(nB,1);
        indInValue_rvs(bamp_in_reg) = indInUni_rvs(bactData.indInSeqs(bamp_in_reg,rr));
        Ad_rvs = sparse(indInValue_rvs(bamp_in_reg),bamp_in_reg,ones(1,length(bamp_in_reg)),Kd,nB);
        
        % Write the matrix
        dat0.M{rr} = [Ad_fwd;Ad_rvs];
        dat0.kmers{rr} = [values_fwd; values_rvs];
    end
    sumMPerBactPerRegion = sum(dat0.M{rr},1);
    dat0.M{rr} = bsxfun(@rdivide,dat0.M{rr},(sumMPerBactPerRegion+eps));
    
    
    % Count errors between reads and kmers
    if AlgoConfig.verbose
        disp('Building matrix Q')
    end
    nY = length(reads(rr).uniqueReads_count);
    cell_A_i = cell(1,nY);
    cell_A_j = cell(1,nY);
    cell_A_s = cell(1,nY);
    non_zero_counts = zeros(1,nY);

    nK = size(dat0.M{rr},1);
    debug_dist_mat = zeros(nY,size(dat0.kmers{rr},1));
    rrL = size(dat0.kmers{rr},2);
    for yy = 1:nY
        if AlgoConfig.verbose && yy == 1000*fix(yy/1000)
            disp(['Mapped ' num2str(yy) ' reads out of ' num2str(nY)])
        end
        
        distvec = sum(bsxfun(@ne,dat0.kmers{rr},reads(rr).uniqueReads(yy,:)),2)-sum(reads(rr).uniqueReads(yy,:)=='N');
        debug_dist_mat(yy,:) = distvec;
        min_dist_2db = min(distvec);
            
        if min_dist_2db <= nMM_cut
            
            mapped_kmer = distvec<=nMM_cut;
            tmp_Pr_read_given_kmer = ((1-pe).^(rrL-distvec(mapped_kmer))).*((pe/3).^distvec(mapped_kmer));
            Pr_read_given_kmer = sparse(find(mapped_kmer),ones(sum(mapped_kmer),1),tmp_Pr_read_given_kmer,nK,1);
            Pr_read_given_j = Pr_read_given_kmer'*dat0.M{rr};
            
            mapped_bact = Pr_read_given_j>0;
            cell_A_i{yy} = yy*ones(1,sum(mapped_bact));
            cell_A_j{yy} = find(mapped_bact);
            cell_A_s{yy} = full(Pr_read_given_j(mapped_bact));
            non_zero_counts(yy) = sum(mapped_bact);
        end
        
    end
    
    
    % Build one sparse matrix
    i_all = [cell_A_i{:}];
    j_all = [cell_A_j{:}];
    s_all = [cell_A_s{:}];
    dat0.A{rr} = sparse(i_all,j_all,s_all,nY,nB);
    dat0.reads{rr} = reads(rr).uniqueReads;
    dat0.F{rr} = reads(rr).uniqueReads_count;
    
    disp('--------------------------------------------')
end
dat0.is_perfect_match = bactData.is_perfect_match;
dat0.indInSeqs = bactData.indInSeqs;

