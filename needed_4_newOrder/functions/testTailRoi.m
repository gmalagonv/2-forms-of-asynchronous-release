function maskSyn = testTailRoi (syn, totNumFrames, movie, lastFrameBeforeTail)



fstFramePart = lastFrameBeforeTail -10;
framesPart = (fstFramePart : 1 : totNumFrames);
finalLengthpart =  length(framesPart);

% [concatpresent] = checkIfConcInFile(nameMov);
% if concatpresent
% %[volume, numImages] = bfOpen3DVolume_GM(filename, avgFlag, frames2load, quietFlag,varargin)
%     [movie, ~ ]= bfOpen3DVolume_GM(nameMov, 0, framesPart,0);
% else
%     %disp('file NO concatenated, load by parts')
%     [movie,~] = movieMergerMTLB(nameMov, framesPart);
% end

movie = movie(:,:,(framesPart));

maskSyn = true(length(syn),1);

for  i = 1:length(syn)
    pieceMovie = reshape(movie((syn(i).rowNumCorr -20 : syn(i).rowNumCorr + 20),(syn(i).columnNumCorr - 20 : syn(i).columnNumCorr + 20),:),[41,41, finalLengthpart]);
    temp = double(syn(i).mask).*double(pieceMovie);
    area = sum(sum(double(syn(i).mask)));
    pieceTrace = squeeze(double(sum(sum(temp))/area));
    
    f0 = mean(pieceTrace(1:10));% last 10 points before tail
    SDtrace = std(pieceTrace(1:10));
    tail = pieceTrace(11:length(pieceTrace));
    tail = denoiseForwrdBack_g(tail,20,10);
    
    [~, locMax] = max(tail);
    if locMax < 4
        locMax = 4;
    end

    if (locMax + 3) > length(tail)
        locMax = length(tail) - 3;
    end
    f = mean(tail((locMax - 3):(locMax + 3))); %% 6 pnts
    %df = (f-f0)/f0;
    maskSyn(i) = (f > (f0 + 2*SDtrace));

    
    
end 



end