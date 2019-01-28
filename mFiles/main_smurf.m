function main_smurf(config_file)

% *********************** LOAD METHOD's PARAMETERS *******************
if nargin == 0
    if ispc
        [ConfigFile, ConfigPath] = uigetfile({'*.m','Matlab Files';'*.txt','Text Files'},'Select Config File');
        if ConfigFile == 0
            return
        end
        config_file = fullfile(ConfigPath,ConfigFile);
    else
        disp('Config file name should be provided as an argument')
        return
    end
end
if ~exist(config_file,'file')
    disp(['Couldn''t find config file: ' config_file])
    return
end
    
% run(fullfile(ConfigPath,ConfigFile));
fid = fopen(config_file);
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    eval(tline)
end
fclose(fid);

% *************************** PROFILE ONE SAMPLE *********************
SampleConfig = struct;
SampleConfig.sample_name = sample_name;
SampleConfig.sample_dir = [base_samples_dir '/' sample_name];
SampleConfig.primers_seq = primers_seq;

main_multiple_regions(PrepConfig,AlgoConfig,SampleConfig)
