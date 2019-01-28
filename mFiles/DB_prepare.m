clear

% Load the primers and pathways
run('../Configs/db_params_script')

packed_sequence_uni = [];
packed_header_uni = [];
len_uni = [];

InfoStruct = fastainfo([db_dir db_filename]);
N = double(InfoStruct.NumberOfEntries);
L = ceil(2*max_s16_len/pack_word_size);
for ii = 1:ceil(N/B)
    disp(['Block '  num2str(ii) '/' num2str(ceil(N/B))])
    [Heads, Seqs] = fastaread([db_dir db_filename],'blockread', [B*(ii-1)+1 B*ii]);
    seq_len = cellfun(@length,Seqs);
    Heads(seq_len>max_s16_len | seq_len<min_s16_len) = [];
    Seqs(seq_len>max_s16_len | seq_len<min_s16_len) = [];
    disp('Done loading. Working...')
    
    F = 2^maxNonACGT;
    iF = 0.1^ceil(log10(4^maxNonACGT+1));
    ind2orig = zeros(1,F*B); sN = 0;   delNum = 0; multNum = 0; new_heads = cell(1,F*B); new_seqs = zeros(F*B, L,['uint' num2str(pack_word_size)]);    new_lens = zeros(1,F*B);
    for ss = 1:length(Seqs)
        if ss == 10000*fix(ss/10000)
            disp(['Counting: ' num2str(ss)])
        end
        [add_seqs,add_heads] = unambiguit_one_seq(Seqs{ss},Heads{ss},maxNonACGT);
        if size(add_seqs,1)>0
            [packed_seqs, packed_lens] = pack_seqs(add_seqs,pack_word_size);
            [pN,pM] = size(packed_seqs);
            
            new_seqs(sN+1:sN+pN,1:pM) = packed_seqs;
            new_lens(sN+1:sN+pN) = packed_lens;
            new_heads(sN+1:sN+pN) = num2cell(str2num(Heads{ss})+iF*(0:pN-1));
            sN = sN + pN;
            
            if size(add_seqs,1)>1
                multNum = multNum + 1;
            end
        else
            delNum = delNum + 1;
        end
    end
    
    % Sort sequences
    [Seqs_sort,xi] = numerical_sort_rows([packed_sequence_uni; new_seqs(1:sN,:)]);
    tmp_Heads = [packed_header_uni new_heads(1:sN)];
    Heads_sort = tmp_Heads(xi);
    tmp_Lens = [len_uni new_lens(1:sN)];
    Lens_sort = tmp_Lens(xi);
    
    % Find unique reads and their headers
    ui = find(sum(abs(diff(double(Seqs_sort),1,1)),2));
    ia = [1; ui+1];
    ib = [ui;size(Seqs_sort,1)];
    wi = find((ia-ib)~=0);
    for ww = 1:length(wi)
        Heads_sort{ia(wi(ww))} = unique([Heads_sort{ia(wi(ww)):ib(wi(ww))}]);
    end
    
    packed_sequence_uni = Seqs_sort(ia,:);
    packed_header_uni = Heads_sort(ia);
    len_uni = Lens_sort(ia);
    
    disp(['Deleted: ' num2str(delNum) ';Ambiguous: ' num2str(multNum) '-became: ' num2str(sN-B-multNum) ';Unique: ' num2str(length(len_uni))])
end
clear Heads Seqs Heads_sort Seqs_sort Lens_sort new_seqs new_heads new_lens

% Dictionary generation
Dic.last_key = max(cellfun(@max,packed_header_uni)); Dic.lists = {}; Dic.headers = {};
[Header_uni, Dictionary] = update_header_dictionary(packed_header_uni,Dic);
Sequence_uni_packed = packed_sequence_uni;


% Prepare the results directory
dot_ind = find(db_filename == '.',1,'last');
suffix = ['unique_up_to_' num2str(maxNonACGT) '_ambiguous_16S'];
dbDirName = [db_dir suffix];
fullname = [db_filename(1:dot_ind-1) '_' suffix];
if exist(dbDirName,'dir')
    rmdir(dbDirName, 's')
end
mkdir(dbDirName)
fa_savename = [dbDirName '/' fullname '.fasta'];
mat_savename = [dbDirName '/' fullname '.mat'];


% Unpack and Save the results to fasta
warning off
N = length(Header_uni);
Sequence_uni = cell(1,N);
Header_uni_str = cell(1,N);
for nn = 1:N
    if nn == 10000*fix(nn/10000)
        disp(['Counting fasta write: ' num2str(nn)])
    end
    
    Sequence_uni{nn} = unpack_seqs(Sequence_uni_packed(nn,:), len_uni(nn), pack_word_size);
    Header_uni_str{nn} = num2str(Header_uni{nn});
    fastawrite(fa_savename, Header_uni_str{nn}, Sequence_uni{nn})
end
warning on


% Save the results to matlab
save(mat_savename,'pack_word_size', 'len_uni', 'Dictionary', 'Header_uni', 'Sequence_uni_packed', 'Sequence_uni', '-v7.3')



