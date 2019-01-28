

% Set the 16S reference DB
uniS16_dir = '../Green_Genes_201305/unique_up_to_3_ambiguous_16S';
db_filename = 'GreenGenes_201305_unique_up_to_3_ambiguous_16S';


% The primers seqs
switch primer_set_name
    case 'amp6Regions'
        % 6 regions primers
        primers_seq{1,1} = 'TGGCGGACGGGTGAGTAA';
        primers_seq{1,2} = 'CTGCTGCCTCCCGTAGGA';
        primers_seq{2,1} = 'TCCTACGGGAGGCAGCAG';
        primers_seq{2,2} = 'TATTACCGCGGCTGCTGG';
        primers_seq{3,1} = 'CAGCAGCCGCGGTAATAC';
        primers_seq{3,2} = 'CGCATTTCACCGCTACAC';
        primers_seq{4,1} = 'AGGATTAGATACCCTGGT';
        primers_seq{4,2} = 'GAATTAAACCACATGCTC';
        primers_seq{5,1} = 'GCACAAGCGGTGGAGCAT';
        primers_seq{5,2} = 'CGCTCGTTGCGGGACTTA';
        primers_seq{6,1} = 'AGGAAGGTGGGGATGACG';
        primers_seq{6,2} = 'CCCGGGAACGTATTCACC';

        DB_kmer_len = 130;
    case 'V3V4'
        % V3-V4 regions primers
        primers_seq{1,1} = 'CCTACGGGNGGCWGCAG';
        primers_seq{1,2} = 'GACTACHVGGGTATCTAATCC';

        DB_kmer_len = 240;
    case 'V4'
        % V4 regions primers
        primers_seq{1,1} = 'GTGCCAGCMGCCGCGGTAA';
        primers_seq{1,2} = 'GGACTACHVGGGTWTCTAAT';

        DB_kmer_len = 240;
end

% Generate kmers DB path and filename
allowed_mm = 2;
suffix = ['_' primer_set_name '_' num2str(allowed_mm) 'mm_RL' num2str(DB_kmer_len)];


