 function[localBckgMov, trace] = localBck(DurationGroupsBckg,movie,InfoExp, varargin)
 %% function 2 define local bckg as the MEDIAN of the prestim frames in a given time interval
%DurationGroupsBckg = 5; %% needed when localBckgFlag == 1 
factorCorr = 0.95;

framesAftStim = length(InfoExp.descartFrame);
numBckg = ceil(InfoExp.NumberAP/DurationGroupsBckg);
[~, firstPreframe] = max(InfoExp.descartFrame);
%mean_noresp_dt = zeros(41,41,numBckg);
localBckgMov = zeros(41,41,length(movie));
fr = 1;
%indxLocalBckg =  ceil(frame/(DurationGroupsBckg*20));

for gr = 1:numBckg
    
    factor = firstPreframe + (gr - 1) * DurationGroupsBckg * framesAftStim;
    localBckgIndx = (factor:framesAftStim:(factor + (DurationGroupsBckg - 1) * framesAftStim));
    %localBckgMov(:,:,(fr:localBckgIndx(end)+19)) = median(movie(:, :, localBckgIndx),3);
    %fr = localBckgIndx(end) + 20;
    if gr ~= numBckg
        for i = fr: (localBckgIndx(end) + framesAftStim -1)
            localBckgMov(:,:,i) = (median(movie(:, :, localBckgIndx),3)).* factorCorr;
        end
    else
         for i = fr: length(movie)
             localBckgMov(:,:,i) = median(movie(:, :, localBckgIndx),3).* factorCorr;
         end
    end
    fr = localBckgIndx(end) + framesAftStim;

end
if ~isempty (varargin) 
    mask = varargin{1};
    temp = (mask).* double(localBckgMov);
    area =  sum(sum(mask));
    trace = squeeze(double(sum(sum(temp))/area));
else
    trace = 0;
end

localBckgMov = cast(localBckgMov, class(movie));
end
