function read_fastq_save_unireads(fname, PrepConfig, dest_pe_flag, readsStatsObj)

% Extract params
data_type = PrepConfig.data_type;
qual_th = PrepConfig.qual_th;
prc_high_qual_th = PrepConfig.prc_high_qual;
low10_th = PrepConfig.low10_th; 
read_len = PrepConfig.read_len; 
pe_flag = PrepConfig.pe_flag; 
max_num_Ns = PrepConfig.max_num_Ns;


%     the 25 prcentile should be hogher than 30, and number of bases
%     smaller than 10 is less than 3
B = 1e5;
Suni = [];
freq = [];

slsh_ind = find(fname == '/',1,'last');
tmp_names = dir([fname '*.' data_type]);

% Dealing with compression due to data2 storage problems
gz_names = dir([fname '*.fastq.gz']);
if isempty(tmp_names) && ~isempty(gz_names)
    for gz = 1:length(gz_names)
        gunzip([fname(1:slsh_ind) gz_names(gz).name])
    end
    tmp_names = dir([fname '*.' data_type]);
end

if pe_flag == 1
    % Match paired end files
    pairs_ind = zeros(0,2);
    for f1 = 1:length(tmp_names)
        tmp_str = tmp_names(f1).name;
        r1_ind = strfind(tmp_str,'_R1_');
        if isempty(r1_ind)
            continue
        else
            tmp_str(r1_ind+2) = '2';
            match = 0;
            for f2 = 1:length(tmp_names)
                if strcmp(tmp_names(f2).name,tmp_str)
                    match = 1;
                    break;
                end
            end
            if match == 1
                pairs_ind(end+1,:) = [f1 f2];
            else
                error('Couldnt match R1-R2 pair')
            end
        end
    end
else
    pairs_ind = (1:length(tmp_names))';
end

% Read the files in pairs and generate unique list of reads
nReads = 0;
nLongReads = 0;
nGoodReads = 0;
for pp = 1:size(pairs_ind,1)
    fname_r1 = [fname(1:slsh_ind) tmp_names(pairs_ind(pp,1)).name];
    if pe_flag == 1
        fname_r2 = [fname(1:slsh_ind) tmp_names(pairs_ind(pp,2)).name];
    end
    
    if strcmp(data_type,'fastq')
        InfoStruct = fastqinfo(fname_r1);
    elseif strcmp(data_type,'fasta')
        InfoStruct = fastainfo(fname_r1);
    end
    N = double(InfoStruct.NumberOfEntries);
    for ii=1:ceil(N/B)
        disp(['Part ' num2str(pp) '/' num2str(size(pairs_ind,1)) ' - Block '  num2str(ii) '/' num2str(ceil(N/B))])
        
        % Filter reads less than read_len and cut reads longer than
        % read_len so all reads are the same length
        if strcmp(data_type,'fastq')
            [H1, S1_cell, q1_cell] = fastqread(fname_r1,'blockread', [B*(ii-1)+1 B*ii]);
        elseif strcmp(data_type,'fasta')
            [H1, S1_cell] = fastaread(fname_r1,'blockread', [B*(ii-1)+1 B*ii]);
        end
        nB = length(S1_cell);

        S1_use_ind = cellfun(@length,S1_cell)>=read_len;
        
        S1 = char(double('N')*ones(nB,read_len));
        S1(S1_use_ind,:) = cell2mat(cellfun(@(x) x(1:read_len),S1_cell(S1_use_ind),'UniformOutput',false)');
        nonambig_1 = sum(S1 ~= 'A' & S1 ~= 'C' & S1 ~= 'G' & S1 ~= 'T',2) <= max_num_Ns;
        
        if strcmp(data_type,'fastq')
            q1 = zeros(nB,read_len);
            q1(S1_use_ind,:) = abs(cell2mat(cellfun(@(x) x(1:read_len),q1_cell(S1_use_ind),'UniformOutput',false)'))-33;
            %         goodqual_1 = prctile(q1,25,2)>prc25_th & sum(q1<10,2)<low10_th;
            goodqual_1 =  mean(q1>qual_th,2)>prc_high_qual_th-eps & sum(q1<10,2)<low10_th;
        else
            goodqual_1 =  true(size(S1,1),1);
        end

        
        % Add this block reads to the list
        if pe_flag == 1
            if strcmp(data_type,'fastq')
                [H2, S2_cell, q2_cell] = fastqread(fname_r2,'blockread', [B*(ii-1)+1 B*ii]);
            elseif strcmp(data_type,'fasta')
                [H2, S2_cell] = fastaread(fname_r2,'blockread', [B*(ii-1)+1 B*ii]);
            end
            S2_use_ind = cellfun(@length,S2_cell)>=read_len;
            
            S2 = char(double('N')*ones(nB,read_len));
            S2(S2_use_ind,:) = cell2mat(cellfun(@(x) x(1:read_len),S2_cell(S2_use_ind),'UniformOutput',false)');
            nonambig_2 = sum(S2 ~= 'A' & S2 ~= 'C' & S2 ~= 'G' & S2 ~= 'T',2) <= max_num_Ns;
            
            if strcmp(data_type,'fastq')
                q2 = zeros(nB,read_len);
                q2(S2_use_ind,:) = abs(cell2mat(cellfun(@(x) x(1:read_len),q2_cell(S2_use_ind),'UniformOutput',false)'))-33;
                %         goodqual_2 = prctile(q2,25,2)>prc25_th & sum(q2<10,2)<low10_th;
                goodqual_2 =  mean(q2>qual_th,2)>prc_high_qual_th-eps & sum(q1<10,2)<low10_th;
            else
                goodqual_2 =  true(size(S2,1),1);
            end

            % Filter bad reads
            S1 = S1(goodqual_1 & goodqual_2 & nonambig_1 & nonambig_2,:);
            S2 = S2(goodqual_1 & goodqual_2 & nonambig_1 & nonambig_2,:);
            if dest_pe_flag == 1
                %     connect pairs
                S = [Suni;[S1,S2]];
                freq = [freq;ones(size(S1,1),1)];
            else
                % use each pair element as single read
                S = [Suni;[S1;S2]];                
                freq = [freq;ones(2*size(S1,1),1)];
            end
            nReads = nReads + nB;
            nLongReads = nLongReads + sum(S1_use_ind & S2_use_ind);
            nGoodReads = nGoodReads + size(S1,1);
        else
            % single end input ==> single end output
            S1 = S1(goodqual_1 & nonambig_1,:);
            S = [Suni;S1];
            freq = [freq;ones(size(S1,1),1)];
            nReads = nReads + nB;
            nLongReads = nLongReads + sum(S1_use_ind);
            nGoodReads = nGoodReads + size(S1,1);
        end            
        
        % Find unique reads and their count
        [S,xi] = sortrows(S);
        freq1 = freq(xi);
        [Suni, ia] = unique(S,'rows', 'first');
        [~, ib] = unique(S,'rows', 'last');
        cumcount = [0;cumsum(freq1)];
        freq = cumcount(ib+1)-cumcount(ia);

    end
 
    % Add stats  
    addStats(readsStatsObj,'Number of loaded reads', NaN, nReads);
    disp(['Number of reads: ' num2str(nReads)])
    
    addStats(readsStatsObj,'Number of long reads', NaN, nLongReads);
    disp(['Percent of long enough reads: ' num2str(nLongReads/nReads)])
    
    addStats(readsStatsObj,'Number of good reads', length(freq), nGoodReads);
    disp(['Percent of good reads: ' num2str(nGoodReads/nLongReads)])
    
end


savename_1 = [fname,'_R1_unireads.fasta'];
if exist(savename_1,'file')
    delete(savename_1)
end
if dest_pe_flag == 1
    savename_2 = [fname,'_R2_unireads.fasta'];
    if exist(savename_2,'file')
        delete(savename_2)
    end
end


tic
B = 100000;
N = size(Suni,1);
n_dgt = ceil(log10(N));
Huni = num2str((1:N)',sprintf('%%%d.%di',n_dgt,n_dgt));
for nn = 1:ceil(N/B)
    disp(['Counting fasta write: ' num2str(nn)])
    fastawrite(savename_1, Huni((nn-1)*B+1:min(nn*B,N),:), Suni((nn-1)*B+1:min(nn*B,N),1:read_len))
    if dest_pe_flag == 1
        fastawrite(savename_2, Huni((nn-1)*B+1:min(nn*B,N),:), Suni((nn-1)*B+1:min(nn*B,N),read_len+1:end))
    end
    toc
end

% Save the unique reads 
save([fname,'_unireads.mat'],'freq','Suni','Huni')
