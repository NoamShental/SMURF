function save_reconstruction(readsPath, sample_num, rdp_db_dir)

if ~exist([rdp_db_dir '/descriptions'],'dir') || ~exist([rdp_db_dir '/key.mat'],'file')
    error('Wrong RDP descriptions directory')
end

% Load the algorithm results
matlab_filename = [readsPath '/resDir/sample_' sample_num '_results.mat'];
load(matlab_filename, 'found_bacteria','bactMetaGroups')



% **************** Build the RECONSTRUCTIONS file ************
[Groups,levels_list] = read_groups_descriptions(rdp_db_dir, found_bacteria, bactMetaGroups);

output_cell = [{'Group index','Frequency', 'Read count', '# of sequences in group'}, levels_list, {'Fraction of RDP matches'}];
output_cell = [output_cell;cat(1, Groups.output)];

% Save to file
reconst_filename = [readsPath '/resDir/sample_' sample_num '_reconstruction.txt'];
saveCellFile(output_cell, reconst_filename)
% ----------------------------------------------------------------------------------



% **************** Build the MAT file ************
groups_filename = [readsPath '/resDir/sample_' sample_num '_reconstruction.mat'];
save(groups_filename,'Groups','output_cell','levels_list')
% ----------------------------------------------------------------------------------

