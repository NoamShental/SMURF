function [input_mat, I] = numerical_sort_rows(input_mat)

[sorted_first_col, I] = sort(input_mat(:,1));
input_mat = input_mat(I,:);

if size(input_mat,2) > 1
    diff_ind = find(diff(sorted_first_col));
    work_ind = [[1;diff_ind+1] [diff_ind;length(sorted_first_col)]];
    work_ind((work_ind(:,1)-work_ind(:,2))==0,:) = []; 
    
    for ii = 1:size(work_ind,1)
        curr_ind = (work_ind(ii,1):work_ind(ii,2))';
        [tmp_out,tmp_I] = numerical_sort_rows(input_mat(curr_ind,2:end));
        input_mat(curr_ind,2:end) = tmp_out;
        I(curr_ind) = I(curr_ind(tmp_I));
    end
end

