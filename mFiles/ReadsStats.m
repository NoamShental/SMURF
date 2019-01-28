classdef ReadsStats < dynamicprops
    properties
        nR 
        
        row_headers = {}
        stats_reads
        stats_counts
        
        regions_row_headers = {}
        regions_stats_reads
        regions_stats_counts
    end
    methods
        function [obj] = ReadsStats(nR)
            obj.nR = nR;
        end
        function addStats(obj,desc, reads_stats, count_stats, rr)
            if ~exist('rr','var')
                rr = 0;
            end
            if length(reads_stats) == 1
                if rr == 0
                    obj.row_headers{end+1,1} = desc;
                    obj.stats_reads = [obj.stats_reads;reads_stats];
                    obj.stats_counts = [obj.stats_counts;count_stats];
                elseif rr<=obj.nR
                    [~,loc] = ismember(desc,obj.regions_row_headers);
                    if loc == 0
                        obj.regions_row_headers{end+1,1} = desc;
                        loc = length(obj.regions_row_headers);
                    end
                    obj.regions_stats_reads(loc,rr) = reads_stats;
                    obj.regions_stats_counts(loc,rr) = count_stats;
                else
                    error('Region number is larger than number of regions')
                end
            elseif length(reads_stats) == obj.nR
                obj.row_regions_headers{end+1,1} = desc;
                obj.regions_stats_reads = [obj.regions_stats_reads;reads_stats];
                obj.regions_stats_counts = [obj.regions_stats_counts;count_stats];
            else
               error('Invalid stasts dimension')
            end
        end
        function disp(obj)
            for st = 1:length(obj.row_headers)
                disp([obj.row_headers{st} ' = ' num2str(obj.stats_counts(st))])
            end
            for st = 1:length(obj.regions_row_headers)
                disp([obj.regions_row_headers{st} ': ' num2str(obj.regions_stats_counts(st,:))])
            end
        end
        function [stat_value] = getOneStat(obj,stat_name)
            [~,loc] = ismember(stat_name,obj.row_headers);
            if loc == 0
                error('Stat name not found')
            end
            stat_value = obj.stats_counts(loc);
        end
    end
end
