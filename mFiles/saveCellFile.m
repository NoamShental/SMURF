%This function saves cell table into outFile 
%function saveCellFile(table, outFile)
function saveCellFile(table, outFile, delimiter)

if ~exist('delimiter','var') || isempty(delimiter)
    delimiter = '\t';
end

fout=fopen(outFile, 'wt');

for i=1:size(table, 1)
    for j=1:size(table, 2)-1
        if isempty(table{i,j})
            fprintf(fout, delimiter);
            continue;
        elseif isnumeric(table{i,j})
            specifier=['%f' delimiter];
        else
            specifier=['%s' delimiter];
        end
        fprintf(fout, specifier, table{i,j});
    end
    if isnumeric(table{i, size(table, 2)})
        specifier='%f';
    else
        specifier='%s';
    end
    fprintf(fout, specifier, table{i, size(table, 2)} );
    fprintf(fout, '\n');
end
fclose(fout);