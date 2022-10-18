function [output] = checkIfConcInFile (nameF)
nameF = strsplit(nameF,'.');
name = nameF{1};
ext = ['.' nameF{2}];

listing = dir;
ALLnames = ({listing(:).name})';
EXTidx = strfind(ALLnames, ext);
EXTidx = cellfun(@(x) ~isempty(x),EXTidx,'UniformOutput',false);
EXTnames = ALLnames(cell2mat(EXTidx));
EXTnames = cellfun(@(x)regexprep(x, ext, ''),EXTnames,'UniformOutput',false);

nameIndx = strfind(EXTnames, name);

%nameIndx = strfind(EXTnames, '_X');
nameIndx = cellfun(@(x) ~isempty(x), nameIndx,'UniformOutput',false);
output = logical(sum(cell2mat(nameIndx)));



