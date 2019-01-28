function save_reconstruction_github(readsPath, sample_name, Taxonomy, Header_uni, Sequence_uni)


% Load the algorithm results
matlab_filename = [readsPath '/resDir/sample_' sample_name '_results.mat'];
load(matlab_filename, 'found_bacteria','bactMetaGroups')

output_cell = {'Group Index', 'Frequency', 'Read count','# of seqs in group'};
if ~isempty(Taxonomy)
    L = length(Taxonomy.ranks_to_extract);
    output_cell = [output_cell {'Fraction of sequences'} Taxonomy.ranks_to_extract];
else
    L = 0;
end

% **************** Build the RECONSTRUCTION groups **************
Groups = struct('freq', [], 'reads', [], 'fractions', [], 'hash', [], 'headers', [], 'answer_cell', []);
Groups(1) = [];
for uu = 1:length(found_bacteria.frequency)
    
    Groups(end+1).freq = found_bacteria.frequency(uu);
    Groups(end).reads = found_bacteria.assigned_reads(uu);
    bact_ind = sort(bactMetaGroups(uu).db_ind);
    
    if L > 0
        this_group_taxa_table = cell2table(Taxonomy.taxa_name_calls(bact_ind,:));
        for ll = 1:L
            [U,~,J] = unique(this_group_taxa_table(:,1:ll),'rows');
            full_names = table2cell(U);
            N = hist(J,1:size(U,1));
            fractions = N'/sum(N);
            Groups(end).fractions{ll} = fractions;
            Groups(end).answer_cell{ll} = [cellstr(num2str(Groups(end).freq*fractions)),cellstr(num2str(round(Groups(end).reads*fractions))), full_names];
        end
        nrows = length(fractions);
        reconstruction_cell = [num2cell([uu*ones(nrows,1), Groups(end).freq*ones(nrows,1), Groups(end).reads*ones(nrows,1), length(bact_ind)*ones(nrows,1), fractions]) full_names];
    else
        reconstruction_cell = num2cell([uu, Groups(end).freq, Groups(end).reads, length(bact_ind)]);
    end
    
    output_cell = [output_cell; reconstruction_cell];
end


% % **************** Save the ANSWER file ************
groups_filename = [readsPath '/resDir/sample_' sample_name '_reconstruction_github.txt'];
saveCellFile(output_cell, groups_filename)


% ******************* Save the MAT file   *************
groups_filename = [readsPath '/resDir/sample_' sample_name '_reconstruction_github.mat'];
save(groups_filename,'Groups')
% ----------------------------------------------------------------------------------

if L == 0
    % **************** Write GROUPS sequences to disk ************
    if ~exist([readsPath '/resDir/groups'],'dir')
        mkdir([readsPath '/resDir/groups'])
    else
        delete([readsPath '/resDir/groups/*.fasta'])
    end
    
    % Write the sequences fasta (group by group)
    if ~isempty(bactMetaGroups)
        for ii = 1:length(bactMetaGroups)
            fasta_filename = [readsPath '/resDir/groups/sample_' sample_name '_group' num2str(ii) '.fasta'];
            fastawrite(fasta_filename, cellfun(@num2str,Header_uni(bactMetaGroups(ii).db_ind),'UniformOutput',false),Sequence_uni(bactMetaGroups(ii).db_ind))
        end
    end
    % ----------------------------------------------------------------------------------
end

