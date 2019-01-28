function main_multiple_regions(PrepConfig,AlgoConfig,SampleConfig,Header_uni,Sequence_uni,Taxonomy)

% Configs processing
sample_name = SampleConfig.sample_name;
readsPath = SampleConfig.sample_dir;
primers_seq = SampleConfig.primers_seq;


% Add primers handling params to AlgoConfig
AlgoConfig.primers_len = PrepConfig.with_primer_flag*zeros(size(primers_seq)) + (1-PrepConfig.with_primer_flag)*cellfun(@length,primers_seq);
AlgoConfig.with_primer_flag = PrepConfig.with_primer_flag;
if ~isfield(AlgoConfig,'use_regions')
    AlgoConfig.use_regions = (1:size(primers_seq,1));
end


% Reads stats struct
readsStatsObj = ReadsStats(size(primers_seq,1));


% Sample preparation
reads_preparation_func([readsPath '/' sample_name], sample_name, primers_seq, PrepConfig, readsStatsObj);

% Save the read count statistics
matlab_filename = [readsPath '/sample_' sample_name '_sampPrepReadStats.mat'];
save(matlab_filename, 'readsStatsObj')


% Community reconstruction
reconstruction_func(readsPath, sample_name, AlgoConfig, readsStatsObj);


% Save the reconstruction
save_reconstruction_github(readsPath, sample_name, Taxonomy, Header_uni, Sequence_uni);

% Save the read count statistics
matlab_filename = [readsPath '/resDir/sample_' sample_name '_readStats.mat'];
save(matlab_filename, 'readsStatsObj')


% Save the reconstruction parameters
matlab_filename = [readsPath '/sample_' sample_name '_paramsStructs.mat'];
save(matlab_filename, 'PrepConfig','AlgoConfig','SampleConfig')
