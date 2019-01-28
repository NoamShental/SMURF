function save_Amat_one_region(filename, primers_seq, Sequence_uni, Header_uni, Dictionary, read_len, allowed_mm, with_primer)

if ~exist('with_primer','var')
    with_primer = 1;
end

rr = str2num(filename(end-4));
nB = length(Sequence_uni);

% -------------------- Forward primer ----------------------
tS = tic;
is_perfect_match_fwd = zeros(nB,size(primers_seq{1},1));
match_pos_fwd = nan(nB,size(primers_seq{1},1));
for aa = 1:size(primers_seq{1},1)
    fp = primers_seq{1}(aa,:); flen = length(fp);
    
    for ss = 1:nB
        if any(is_perfect_match_fwd(ss,:))
            continue % if already found
        end
        if ss == 1e5*fix(ss/1e5)
            disp(['Region: ' num2str(rr) ' forward'])
            disp(['Primer: ' num2str(aa) ' out of ' num2str(size(primers_seq{1},1))])
            disp(['Amplified: ' num2str(ss) ' sequences out of ' num2str(nB) ])
            toc(tS)
            disp('-----------------------------------')
        end
        
        % Fwd primer match
        Start_fwd = [];
        kfp = strfind(Sequence_uni{ss}, fp);
        if isempty(kfp) && allowed_mm > 0
            [~, Alignment, Start_fwd] = swalign(Sequence_uni{ss}, fp, 'ALPHABET','NT');
            if sum(Alignment(2,:)=='|')>=flen-allowed_mm && sum(sum(Alignment([1 3],:)=='-')) == 0
                kfp = Start_fwd(1) - Start_fwd(2) + 1;
            end
        else
            is_perfect_match_fwd(ss,aa) = 1;
        end
        if ~isempty(kfp)
            match_pos_fwd(ss,aa) = kfp(1);
        end
    end
end

% -------------------- Reverse primer ----------------------
is_perfect_match_rvs = zeros(nB,size(primers_seq{2},1));
match_pos_rvs = nan(nB,size(primers_seq{2},1));
for bb = 1:size(primers_seq{2},1)
    rp = seqrcomplement(primers_seq{2}(bb,:)); rlen = length(rp);
    
    for ss = 1:nB
        if any(is_perfect_match_rvs(ss,:)) 
            continue % if already found
        end
        if ss == 1e5*fix(ss/1e5)
            disp(['Region: ' num2str(rr) ' reverse'])
            disp(['Primer: ' num2str(bb) ' out of ' num2str(size(primers_seq{2},1)) ])
            disp(['Amplified: ' num2str(ss) ' sequences out of ' num2str(nB) ])
            toc(tS)
            disp('-----------------------------------')
        end
        
        % Rvs primer match
        Start_rvs = [];
        krp = strfind(Sequence_uni{ss}, rp);
        if isempty(krp) && allowed_mm > 0
            [~, Alignment, Start_rvs] = swalign(Sequence_uni{ss}, rp, 'ALPHABET','NT');
            if sum(Alignment(2,:)=='|')>=rlen-allowed_mm && sum(sum(Alignment([1 3],:)=='-')) == 0
                krp = Start_rvs(1) - Start_rvs(2) + 1;
            end
        else
            is_perfect_match_rvs(ss,bb) = 1;
        end
        krp = krp - 1;
        if ~isempty(krp)
            match_pos_rvs(ss,bb) = krp(1);
        end
    end
end


% -------------------- Get the sequences ----------------------
[Sequence_amp_mat(1:nB,1:read_len*2)] = deal('N');
is_perfect_match = any(is_perfect_match_fwd,2) & any(is_perfect_match_rvs,2);
for ss = 1:nB
    if ss == 1e5*fix(ss/1e5)
        disp(['Region: ' num2str(rr) ' sequences'])
        disp(['Amplified: ' num2str(ss) ' sequences out of ' num2str(nB)])
        toc(tS)
        disp('-----------------------------------')
    end
    
    if any(is_perfect_match_fwd(ss,:))
        degen_ind_fwd = find(is_perfect_match_fwd(ss,:));
    else
        degen_ind_fwd = find(match_pos_fwd(ss,:)~=0 & ~isnan(match_pos_fwd(ss,:)),1);
    end
    kfp = match_pos_fwd(ss,degen_ind_fwd);

    if any(is_perfect_match_rvs(ss,:))
        degen_ind_rvs = find(is_perfect_match_rvs(ss,:));
    else
        degen_ind_rvs = find(match_pos_rvs(ss,:)~=0 & ~isnan(match_pos_rvs(ss,:)),1);
    end
    krp = match_pos_rvs(ss,degen_ind_rvs);
    
    if ~isempty(kfp) && ~isempty(krp) && ~isnan(kfp) && ~isnan(krp)
        amp_len = krp(1) - kfp(1) + rlen*with_primer - flen*(1-with_primer) + 1;
        if amp_len >= read_len
            % Fwd read
            if with_primer == 0
                kfp = kfp + flen;
            end
            if kfp>0
                Sequence_amp_mat(ss,1:read_len) = Sequence_uni{ss}(kfp:kfp+read_len-1);
            else
                Sequence_amp_mat(ss,Start_fwd(2):read_len) = Sequence_uni{ss}(Start_fwd(1):Start_fwd(1)+read_len-Start_fwd(2));
            end
            
            % Rvs read
            if with_primer == 1
                krp = krp+rlen;
            end
            sL = length(Sequence_uni{ss});
            if krp<=sL
                Sequence_amp_mat(ss,read_len+1:end) = Sequence_uni{ss}(krp-read_len+1:krp);
            else
                Sequence_amp_mat(ss,read_len+1:end-(krp-sL)) = Sequence_uni{ss}(krp-read_len+1:end);
            end
        end
    end
    
end
Header_amp = cell2mat(Header_uni);

% Find unique reads
[values, iI, indInValue] = unique(Sequence_amp_mat,'rows');

% Remove the not found
NNN_ind = find(sum(values=='N',2)==2*read_len);
if ~isempty(NNN_ind)
    values(NNN_ind,:) = [];
    indInValue(indInValue==NNN_ind) = 0;
    indInValue(indInValue>NNN_ind) = indInValue(indInValue>NNN_ind) - 1;
end


% Save the database
save(filename,'Header_amp','Dictionary','values','indInValue','is_perfect_match','match_pos_fwd','match_pos_rvs')
