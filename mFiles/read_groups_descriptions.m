function [Groups,levels_list] = read_groups_descriptions(rdp_db_dir, found_bacteria, bactGroups)

% Build the levels list
use_levels = {'domain','phylum','class','order','family','genus','seq_names'};
levels_list = use_levels{1};
for ll = 2:length(use_levels)
    levels_list = [levels_list '::' use_levels{ll}];
end


% Load the original key
load([rdp_db_dir '/key.mat'],'col_names','header2file_map')


% Combine the descriptions per group
Groups = struct('group_cell',[], 'output', [], 'freq', [], 'answer_cell', [], 'name', []);
Groups(1) = [];
SortedGroups = Groups;
n_bact = length(found_bacteria.frequency);
for uu = 1:n_bact
    % Combine all sequences of one group
    bact_ind = sort(bactGroups(uu).db_ind);
    disp(['Loading RDP data of group ' num2str(uu) ' out of ' num2str(n_bact) ': ' num2str(length(bact_ind)) ' sequences'])
    
    len_db_ind = length(bact_ind);
    select_ind = sort(randperm(len_db_ind,min(len_db_ind,1000)));
    oneGroup = read_one_group_descriptions(bact_ind(select_ind),header2file_map,rdp_db_dir);
    
    % Generate the reconstruction description for the whole group
    [full_names,fractions] = build_full_names(oneGroup,use_levels);
    
    % Build the groups cell
    group_params = {num2str(uu), num2str(found_bacteria.frequency(uu)), num2str(found_bacteria.assigned_reads(uu)), num2str(length(bact_ind))};
    Groups(end+1).group_cell = [repmat(group_params,length(fractions),1) full_names fractions];
    Groups(end).output = [Groups(end).group_cell; cellstr(char(9*ones(size(Groups(end).group_cell,2),1)))'];
    
    % Assign a name to the group
    if ~isempty(oneGroup)
        [Ug,~,Jg] = unique(cat(1,oneGroup.genus));
        Ng = hist(Jg,1:length(Ug));
        [~,gg] = max(Ng);
        Groups(end).name = Ug(gg);
    else
        Groups(end).name = 'not found';
    end
    
    % Record the frequency
    Groups(end).index = uu;
    Groups(end).freq = found_bacteria.frequency(uu);
    Groups(end).bact_ind = bact_ind;
        
end

% Sort the  based on genus frequency
sort_order = nan(1,length(Groups));
so_ind = 0;
group_names = [Groups.name];
group_freq = [Groups.freq];
[Ur,Ir,Jr] = unique(group_names);
for gr = 1:length(Ir)
    genus_ind = find(Jr==gr);
    genus_freq = group_freq(genus_ind);
    [~,If] = sort(-genus_freq);
    sort_order(so_ind+1:so_ind+length(If)) = genus_ind(If);
    so_ind = so_ind + length(If);
end
Groups = Groups(sort_order);

%*******************************************************************************

