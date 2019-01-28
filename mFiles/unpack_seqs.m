% Unpack sequences such that each 32-bit word contains 16 basepairs. 
% The unpacking also transfers 0-3 to 1-4 !!    
%
% Input: 
% packed_seqs - array of sequences in packed form (2-bit per nucleotide)
% seqs_len - length of each sequence in nucleotides
% word_size - (new!) enable different word sizes (default is 32 bits but supports any even size) 
%
% Output: 
% seqs - sequences in unpacked form (1-4 representing A,C,G,T)
% 
function seqs = unpack_seqs(packed_seqs, seqs_len, word_size, int2nt)

if(~exist('word_size', 'var') || isempty(word_size)) % default is 32 bits
    word_size = 32;
end
half_word_size = word_size/2;

[num_seqs packed_seqs_len] = size(packed_seqs);

% seqs_len in the input must be between .. (packed_seqs_len-1)*16+1 to packed_seqs_len*16;
% seqs_len = packed_seqs_len*16;

% Here we multiply by a 16 factor in space
seqs = zeros(num_seqs, half_word_size*packed_seqs_len);

% Here we loop only on the 16 places - saves time ...
for i=1:half_word_size
    % Int version
    seqs(:, i:half_word_size:i+half_word_size*(packed_seqs_len-1)) = ...
        uint32(bitget(packed_seqs, 2*i-1)) + 2* uint32(bitget(packed_seqs, 2*i))+1;        
%    seqs(:, i:half_word_size:i+half_word_size*(packed_seqs_len-1)) = mod( floor(packed_seqs / 4^(i-1)), 4 ) + 1;
end
    

seqs = seqs(:, 1:seqs_len); % 'Chop' the unneccessary tail 

if nargin < 4 || int2nt
    seqs(seqs==1) = 'A';
    seqs(seqs==2) = 'C';
    seqs(seqs==3) = 'G';
    seqs(seqs==4) = 'T';
    seqs = char(seqs);
end



