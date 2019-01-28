% Pack the sequences such that every 32-bit word contain 16 basepairs.
% The packing also transfers 1-4 to 0-3 !!
% Note : It is assumed that all sequences here are of the same length here
%
% Input:
% seqs - an array of sequences, in numeric format
% word_size - (new!) enable different word sizes (default is 32 bits but supports any even size)
%
% Output:
% packed_seqs - the sequences in a dense format, 2 bits per nucleotide
% seqs_lens - the length of each sequence
%
function [packed_seqs, seqs_len] = pack_seqs(seqs, word_size)

if(~exist('word_size', 'var') || isempty(word_size)) % default is 32 bits
    word_size = 32;
end
half_word_size = word_size/2;

if(iscell(seqs)) % New! get the function working also for cell array
    n = length(seqs);
    packed_seqs = cell(n,1); seqs_len = zeros(n,1);
    for i=1:n
        [packed_seqs{i} seqs_len(i)] = pack_seqs(seqs{i}, word_size);
    end
else
    [num_seqs seqs_len] = size(seqs);
    if(~isnumeric(seqs)) % convert txt nucleotides to numeric
        %         seqs = nt2int(seqs);
        
        %--------- GARI ------------------
        seqs2 = zeros(size(seqs),'uint8');
        seqs2(seqs=='A') = 1;
        seqs2(seqs=='C') = 2;
        seqs2(seqs=='G') = 3;
        seqs2(seqs=='T') = 4;
        
        if sum(sum(seqs2 == 0))
            error('Packing is supported for ACGT only')
        end
        seqs = seqs2;
        %----------------------------------
    end
    
    
    packed_seqs_len = ceil(seqs_len/half_word_size);  % Here we save a 16 factor in space
    if(word_size == 32)
        seqs = uint32(seqs);
        packed_seqs = zeros(num_seqs, packed_seqs_len, 'uint32');
        base_vec = uint32(4) .^ uint32([0:(half_word_size-1)]); % powers of four
    else
        seqs = uint64(seqs); % here they cosume lots of memory ...
        packed_seqs = zeros(num_seqs, packed_seqs_len, 'uint64');
        base_vec = uint64(4) .^ uint64([0:(half_word_size-1)]); % powers of four
    end
    base_vec = repmat(base_vec, num_seqs, 1);
    for i=1:packed_seqs_len-1 % The last block is dealt with seperately ..
        %         tmp_seqs = (seqs(:,half_word_size*(i-1)+1:half_word_size*i)-1) .* base_vec;
        %         packed_seqs(:,i) = tmp_seqs(:,1) + tmp_seqs(:,2);
        packed_seqs(:,i) = ...
            sum( (seqs(:,half_word_size*(i-1)+1:half_word_size*i)-1) .* base_vec, 2, 'native');
    end
    packed_seqs(:,packed_seqs_len) = sum(( seqs(:,half_word_size*(packed_seqs_len-1)+1:seqs_len) - 1) .* ...
        base_vec(:, 1:mod(seqs_len-1, half_word_size)+1), 2, 'native');      % Here we do the last block which might be shorter ...
end


