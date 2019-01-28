function [new_seqs,new_heads] = unambiguit_one_seq(old_seq,old_head,maxNonACGT)

I = find(old_seq ~= 'A' & old_seq ~= 'C' & old_seq ~= 'G' & old_seq ~= 'T');
if length(I) > maxNonACGT
    new_seqs = [];
    new_heads = [];
    return
elseif isempty(I) 
    new_seqs = old_seq;
    new_heads = old_head;
    return
end

amb_letters = old_seq(I);
letters = 'RYKMSWBDHVN';
disamb_strs = {'AG','CT','GT','AC', 'CG','AT','CGT','AGT', 'ACT', 'ACG', 'ACGT'};
lens = cellfun(@length,disamb_strs);

[TF, dic_ind] = ismember(amb_letters,letters);
if sum(TF) ~= length(TF)
    error('Some ambiguous letters are not found in dictionary')
end


% Generate ambiguous letters matrix
amb_rep_mat = zeros(1,0);
head_add_mat = zeros(1,0);
for ii = length(dic_ind):-1:1
    
    % Sequence replacement
    left_mat = repmat(disamb_strs{dic_ind(ii)},size(amb_rep_mat,1),1);
    left_vec = reshape(left_mat,size(amb_rep_mat,1)*lens(dic_ind(ii)),1);
    right_mat = repmat(amb_rep_mat,lens(dic_ind(ii)),1);
    amb_rep_mat = [left_vec right_mat];
    
    % Headers addition
    head_str = ['_' amb_letters(ii) num2str(I(ii), '%4.4i')];
    left_mat = repmat(head_str,size(head_add_mat,1)*lens(dic_ind(ii)),1);
    right_mat = repmat(head_add_mat,lens(dic_ind(ii)),1);
    head_add_mat = [left_mat left_vec right_mat];

end


% Build the new sequences matrix
N = prod(lens(dic_ind));
new_seqs = repmat(old_seq,N,1);
new_seqs(:,I) = amb_rep_mat;
% new_seqs = sm2cell(new_seqs_mat)';

new_heads = [repmat(old_head,N,1) head_add_mat];
% new_heads = sm2cell(new_heads_mat)';
    
