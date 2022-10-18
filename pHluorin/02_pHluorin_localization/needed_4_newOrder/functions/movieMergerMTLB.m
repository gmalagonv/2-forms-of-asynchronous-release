function [concatenated, totNumFrames] = movieMergerMTLB(nameF, varargin)
%% function to load unconcatenated movie. We miust be laready in the directory where movies are
nameF = strsplit(nameF,'.');
%name = nameF{1};
numb = regexp(nameF{1},'\d*','Match');
ext = ['.' nameF{2}];

listing = dir;
ALLnames = ({listing(:).name})';
EXTidx = strfind(ALLnames, ext);
EXTidx = cellfun(@(x) ~isempty(x),EXTidx,'UniformOutput',false);
EXTnames = ALLnames(cell2mat(EXTidx));

if isempty(varargin)
    frames2load = 1;
    %quietFlag = 1;
else
    frames2load = varargin{1};
    %quietFlag = 0;
end

ROOTidx = strfind(EXTnames, ['_' numb{1} ext ]);
ROOTidx = cellfun(@(x) ~isempty(x), ROOTidx,'UniformOutput',false);
EXTnames = cellfun(@(x)regexprep(x, ext, ''),EXTnames,'UniformOutput',false);
ROOTname = EXTnames(cell2mat(ROOTidx));
if isempty(ROOTname)
    disp('check name, movie not found')
    return
    
else
    PARTroot = [ROOTname{1} '_'];
    PARTidx = strfind(EXTnames, PARTroot);
    PARTidx = cellfun(@(x) ~isempty(x),PARTidx,'UniformOutput',false);
    PARTnames = [ROOTname{1}; EXTnames(cell2mat(PARTidx))];
    concatenated = [];
    for j = 1:length(PARTnames)
        %[movie, ~ ]= bfOpen3DVolume_GM(nameMov, 0, framesPart,0);
        %[volume, numImages] = bfOpen3DVolume_GM(filename, avgFlag, frames2load, quietFlag,varargin)
        [part, ~] = bfOpen3DVolume_GM([PARTnames{j} ext], 0, frames2load, 1);
        concatenated = cat(3,concatenated, part);
        
        disp([PARTnames{j} ' loaded'])
        
    end
    totNumFrames = length(concatenated);
end

end


