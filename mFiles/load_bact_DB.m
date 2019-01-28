function [bactData] = load_bact_DB(dbPath,dbFileName,with_primer_flag,primers_len,read_len,const_len_flag,use_regions,prob_amp)

regions_files = dir([dbPath '/' dbFileName '*.mat'])
if ~exist('use_regions','var') || isempty(use_regions)
    use_regions = 1:length(regions_files);
end
nR = length(use_regions);

if ~exist('const_len_flag','var') || isempty(const_len_flag)
    const_len_flag = 0;
end

if ~exist('with_primer_flag','var') || isempty(with_primer_flag) || with_primer_flag == 1
    with_primer_flag = 1;
    primers_len = zeros(nR,2);
end

if ~exist('prob_amp','var') || isempty(prob_amp)
    prob_amp = 1;
end


tmp_ind = strfind(dbFileName,'RL');
kL = num2str(dbFileName(tmp_ind+2:end));
if ~exist('read_len','var') || isempty(read_len)
    read_len = str2num(kL);
end

[dbPath '/' regions_files(1).name]
load([dbPath '/' regions_files(1).name],'Dictionary','Header_amp');
nB = length(Header_amp);
bactData.Header_Dictionary = Dictionary;
bactData.Header_amp = Header_amp;
bactData.kmers = cell(1,nR);
bactData.indInSeqs = zeros(nB,nR);
bactData.is_perfect_match = zeros(nB,nR);

for new_rr = 1:nR
    rr = use_regions(new_rr);
    
    disp(['Loading bacterial DB for region ' num2str(new_rr) ' out of ' num2str(nR) ' from original region ' num2str(rr)])
    load([dbPath '/' regions_files(rr).name],'values','is_perfect_match','indInValue');
    
    
    % Reduce primers
    is_amplified = indInValue>0;
    if const_len_flag == 0  
        kmers_no_primer = [values(:,primers_len(rr,1)+1:read_len) values(:,end-read_len+1:end-primers_len(rr,2))];
    else
        kmers_no_primer = [values(:,primers_len(rr,1)+(1:read_len)) values(:,(end-read_len+1:end)-primers_len(rr,2))];
    end
    [unique_no_primer, ~, indInUni] = unique(kmers_no_primer,'rows');
    indInKmers = zeros(nB,1);
    indInKmers(is_amplified) = indInUni(indInValue(is_amplified));
    
    
    % Reduce non-amplified
    amp_fact = prob_amp/mean(is_amplified);
    new_is_amplified = (rand(nB,1)<amp_fact) & is_amplified;
    
    new_is_perfect_match = is_perfect_match;
    new_is_perfect_match(~new_is_amplified) = false;
    
    new_indInKmers = indInKmers;
    new_indInKmers(~new_is_amplified) = 0;
    
    [keep_kmer_ind,I,J] = unique(new_indInKmers);
    new_unique_no_primer = unique_no_primer(keep_kmer_ind(2:end),:);
    new_indInKmers = J-1;
    
        
    
    bactData.is_perfect_match(:,new_rr) = new_is_perfect_match;
    bactData.kmers{new_rr} = new_unique_no_primer;
    bactData.indInSeqs(:,new_rr) = new_indInKmers;
end

