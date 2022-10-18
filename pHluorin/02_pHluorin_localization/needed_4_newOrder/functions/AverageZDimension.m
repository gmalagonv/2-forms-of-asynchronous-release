function [Average] = AverageZDimension(RawData,starts,ENDs)

if nargin < 3
    
    [~,~,p] = size(RawData); %% p is number of Frames in the movie m x n pixels in the imagin
    Average = sum(RawData,3,'double')/p;
    
else

    N = length(starts);
    RawData2 = [];
    sizeZ = 0;

    for i = 1:N

        start = starts(i);
        END = ENDs(i);
        sample = RawData(:,:,start:END);
        [~, ~, Length] = size(sample);
        RawData2 (:,:,(sizeZ+1):(sizeZ+Length)) = sample;
        [~, ~, sizeZ] = size(RawData);
    end

    [~,~,p] = size(RawData2); %% p is number of Frames in the movie m x n pixels in the imagin
    Average = sum(RawData2,3,'double')/p;
    
   
end
end