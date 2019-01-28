function [names,fractions] = build_full_names(rdp_matches,level_names)

if ~isempty(rdp_matches)
    full_names = cat(1,rdp_matches.(level_names{1}));
    [ddot{1:length(full_names)}] = deal('::');
    ddot = ddot';
    for nn = 2:length(level_names)
        this_level_names = cat(1,rdp_matches.(level_names{nn}));
        TF = ismember(this_level_names,'not found');
        this_level_names(TF) = {['Unknown ' level_names{nn}]};
        
        this_level_names = strcat(ddot,this_level_names);
        full_names = strcat(full_names,this_level_names);
    end
    full_names = cellfun(@(x) strrep(x,'"',''),full_names,'UniformOutput',false);
    
    % Find unique names and their frequencies
    [Uf,~,Jf] = unique(full_names);
    N = hist(Jf,1:length(Uf));
    [~,I] = sort(N,'descend');
    fractions_num = N(I)'/sum(N);
    
    names = Uf(I);
    fractions = cellfun(@num2str,mat2cell(fractions_num,ones(length(I),1),1),'UniformOutput',false');
else
    names = {['Unknown ' level_names{1}]};
    for nn = 2:length(level_names)
        names{1} = [names{1} '::Unknown ' level_names{nn}];
    end
    fractions = {'1'};
end


%*******************************************************************************