function read_unireads_save_split_to_regions(fname,sample_num,primers_seq,max_err_inprimer,with_primer_flag,pe_flag,readsStatsObj)

% broff = 7; % FOR Michael SEcond Lane since the barcode was part of the primer
if ~exist('broff','var')
    broff = 0; % barcode_offset
end
sl_ind = find(fname=='/',1,'last');
sample_dir = fname(1:sl_ind-1);

nR = size(primers_seq,1);
load([fname '_unireads.mat'])
if isempty(Suni)
    frequni = [];
    readsuni = [];
    fwd_rvs_ind = [];
    for rr = 1:nR
        addStats(readsStatsObj,'Number of reads mapped to region', length(frequni), sum(frequni), rr);
        save([sample_dir '/sample_' sample_num '_region_' num2str(rr) '_unireads.mat'],'readsuni','frequni','fwd_rvs_ind');
    end
    return
end

freq_vec = freq;
if pe_flag == 1
    % Reverse the RVS read
    rL = size(Suni,2)/2;
    tmpSuni  = Suni(:,rL+1:end);
    tmpSuniA = Suni(:,rL+1:end);
    tmpSuniA(tmpSuni=='A') = 'T';
    tmpSuniA(tmpSuni=='C') = 'G';
    tmpSuniA(tmpSuni=='G') = 'C';
    tmpSuniA(tmpSuni=='T') = 'A';
    Suni(:,rL+1:end) = fliplr(tmpSuniA);
else
    error('We dont support this I think')
end

mapping_mat = zeros(length(freq_vec),nR);
fwd_mapping_mat = inf(length(freq_vec),nR);
rvs_mapping_mat = inf(length(freq_vec),nR);
for rr = 1:nR
    if pe_flag == 0

        % Forward primers
        for aa = 1:size(primers_seq{rr,1},1)
            flen = size(primers_seq{rr,1},2);
            distvec = sum(bsxfun(@ne,Suni(:,broff+(1:flen)),primers_seq{rr,1}(aa,:)) & Suni(:,broff+(1:flen))~='N',2);
            use_ind = (distvec<=max_err_inprimer);
            mapping_mat(use_ind,rr) = 1;
        end    
        
        % Reverse primers
        for bb = 1:size(primers_seq{rr,2},1)
            rlen = size(primers_seq{rr,2},2);
            distvec = sum(bsxfun(@ne,Suni(:,(end-rlen+1:end)-broff),seqrcomplement(primers_seq{rr,2}(bb,:))) & Suni(:,(end-rlen+1:end)-broff)~='N',2);
            use_ind = (distvec<=max_err_inprimer);
            mapping_mat(use_ind,rr) = -1;
        end
        
        % Get the reads and freqs
        fwd_use_ind = (mapping_mat(:,rr)== 1);
        rvs_use_ind = (mapping_mat(:,rr)==-1);
        if with_primer_flag == 0
            readsuni = [Suni(fwd_use_ind,broff+flen+1:end); Suni(rvs_use_ind,1:end-rlen-broff)];
        else
            readsuni = [Suni(fwd_use_ind,:); Suni(rvs_use_ind,:)];
        end
        frequni = [freq_vec(fwd_use_ind);freq_vec(rvs_use_ind)];
        fwd_rvs_ind = [ones(sum(fwd_use_ind),1); -ones(sum(rvs_use_ind),1)];
        
    elseif pe_flag == 1

        % Mapp the primers
        distvec_fwd = inf(size(Suni,1),1);
        for aa = 1:size(primers_seq{rr,1},1)
            flen = size(primers_seq{rr,1},2);
            distvec_fwd = min(distvec_fwd, sum(bsxfun(@ne,Suni(:,broff+(1:flen)),primers_seq{rr,1}(aa,:)) & Suni(:,broff+(1:flen))~='N' ,2));
        end
        
        distvec_rvs = inf(size(Suni,1),1);
        for bb = 1:size(primers_seq{rr,2},1)
            rlen = size(primers_seq{rr,2},2);
            distvec_rvs = min(distvec_rvs,sum(bsxfun(@ne,Suni(:,(end-rlen+1:end)-broff),seqrcomplement(primers_seq{rr,2}(bb,:))) & Suni(:,(end-rlen+1:end)-broff)~='N',2));
        end
        
        % Assign reads to regions
        use_ind = (distvec_fwd<=max_err_inprimer) & (distvec_rvs<=max_err_inprimer);
        mapping_mat(use_ind,rr) = 1;
        fwd_mapping_mat(:,rr) = distvec_fwd;
        rvs_mapping_mat(:,rr) = distvec_rvs;

        
        % Get the reads and freqs
        pe_use_ind = (mapping_mat(:,rr)==1);
        if with_primer_flag == 0
            readsuni = Suni(pe_use_ind,broff+flen+1:end-rlen-broff);
        else
            readsuni = Suni(pe_use_ind,:);
        end
        frequni = freq_vec(pe_use_ind);
        fwd_rvs_ind = zeros(sum(pe_use_ind),1);

    end
    addStats(readsStatsObj,'Number of reads mapped to region', length(frequni), sum(frequni), rr);
    save([sample_dir '/sample_' sample_num '_region_' num2str(rr) '_unireads.mat'],'readsuni','frequni','fwd_rvs_ind');
end
disp(['Mapped to primers ' num2str(round(100*sum(sum(abs(mapping_mat)))/size(Suni,1))) '% of unique reads'])
disp(['Mapped to primers ' num2str(round(100*sum(freq_vec(sum(abs(mapping_mat),2)>0))/sum(freq_vec))) '% of read counts'])


% [S,xi] = sortrows(Suni(:,end-rlen+1:end));
% freq1 = freq(xi);
% [~, ia] = unique(S,'rows', 'first');
% [~, ib] = unique(S,'rows', 'last');
% cumcount = [0;cumsum(freq1)];
% Pfreq = cumcount(ib+1)-cumcount(ia);
% [Puni, I, J] = unique(S,'rows');

