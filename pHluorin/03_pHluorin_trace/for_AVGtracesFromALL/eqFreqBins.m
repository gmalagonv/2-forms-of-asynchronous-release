function BINS  = eqFreqBins(vector, numBins )
%% function to generate equal frequency bins
% sizeBin = round(length(vector) / numBins);
[~, indxSorted] = sort(vector);
BINS = zeros(size(vector));


    
for i = 1:numBins
    
    BINS((indxSorted(indxSorted(discretize(indxSorted, numBins) == i)))) = i;
%     
%     edges = ((i-1)*sizeBin + 1):(sizeBin*i);
%     if any(edges > length(indxSorted))
%         edges(edges>length(indxSorted)) = length(indxSorted);
%     end
%     indx = indxSorted(edges);
%    BINS(indx) = i;
    
end
