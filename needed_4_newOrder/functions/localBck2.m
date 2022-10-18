 function[localBckgMov, trace] = localBck2(DurationGroupsBckg,medFiltFlag, movie,InfoExp, varargin)
 %% function 2 define local bckg as the MEDIAN of the prestim frames in a given time interval
%DurationGroupsBckg = 5; %% needed when localBckgFlag == 1 

 
if isfield(InfoExp, 'FileExtension') && strcmp(InfoExp.FileExtension,'.tif')
    factorCorr = 0.97;%0.97;
   
else
    factorCorr = 0.96;
end




framesAftStim = length(InfoExp.descartFrame);
numBckg = ceil(InfoExp.NumberAP/DurationGroupsBckg);
[~, firstPreframe] = max(InfoExp.descartFrame);
%mean_noresp_dt = zeros(41,41,numBckg);
localBckgMov = movie;%zeros(41,41,length(movie));
fr = 1;
%indxLocalBckg =  ceil(frame/(DurationGroupsBckg*20));
matrixBckIndx = zeros(framesAftStim,DurationGroupsBckg);
for gr = 1:numBckg
    
    factor = firstPreframe + (gr - 1) * DurationGroupsBckg * framesAftStim;
    startV = factor; 
    
    for j = 1:framesAftStim
        
        matrixBckIndx(j , :) = (startV:framesAftStim:(startV + (DurationGroupsBckg - 1) * framesAftStim));
        startV = startV + 1;
        
    end
    
    localBckgIndx = (factor:framesAftStim:(factor + (DurationGroupsBckg - 1) * framesAftStim));
    %localBckgMov(:,:,(fr:localBckgIndx(end)+19)) = median(movie(:, :, localBckgIndx),3);
    %fr = localBckgIndx(end) + 20;
    if gr ~= numBckg
        
        for k = 1:framesAftStim
            if medFiltFlag == 1
                val = medfilt2(median(movie(:,:,(matrixBckIndx(k,:))),3).*factorCorr);
            else
                val = median(movie(:,:,(matrixBckIndx(k,:))),3).*factorCorr;
            end
 
            for frm = matrixBckIndx(k,:)
                localBckgMov(:,:,frm) = val;
            end

            
        end
        
        
%         for i = fr: (localBckgIndx(end) + framesAftStim -1)
%             localBckgMov(:,:,i) = (median(movie(:, :, localBckgIndx),3)).* factorCorr;
%         end
    else
         for i = fr: length(movie)
             localBckgMov(:,:,i) = median(movie(:, :, localBckgIndx),3).* factorCorr;
         end
    end
    fr = localBckgIndx(end) + framesAftStim;

end

localBckgMov = cat(3, movie(:,:,1), localBckgMov);
localBckgMov(:,:,end) = [];

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
