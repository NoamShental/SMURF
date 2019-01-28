function [oneGroup] = read_one_group_descriptions(use_ind,header2file_map,rdp_db_dir)

oneGroup = [];
for bb = 1:length(use_ind)
    file_index = header2file_map(use_ind(bb),2);
    bb_header = header2file_map(use_ind(bb),1);
    filename = [rdp_db_dir '/descriptions/rdpDescFile' num2str(file_index) '.mat'];
    var_name = ['h' strrep(num2str(bb_header),'.','_')];
    load(filename,var_name)
    if exist(var_name,'var')
        if isempty(oneGroup)
            eval(['oneGroup = ' var_name ';'])
        else
            eval(['oneGroup(end+1) = ' var_name ';'])
        end
        clear(var_name)
    end
end
