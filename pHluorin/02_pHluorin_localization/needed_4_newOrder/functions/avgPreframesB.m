function [preframeMov, trace] = avgPreframesB(movie, meanFlag, randomFlag,varargin)
preframeMov = zeros(size(movie));
factor = 0.9; 
minumPrevFr = 5; 

if randomFlag ==1
    minumPrevFr = 10;
end
    

for fr = 1: length(movie)
    if fr == 1
        preframeMov(:,:,fr) =   movie(:,:,fr);
    else
        if fr <= minumPrevFr
             indexes = (1 : (fr - 1));
        else
            if randomFlag == 1
               indexes = randi([(fr-minumPrevFr) (fr-1)],4,1);
                
                
            else
                indexes = ((fr - minumPrevFr) : 2 :(fr - 1));
            end
        end

       if meanFlag == 1
            preframeMov(:,:,fr) = mean(movie(:, :, indexes),3) .* factor;
        else
            preframeMov(:,:,fr) = median(movie(:, :, indexes),3) .* factor;
        end
            
    end
end

if ~isempty (varargin) 
    mask = varargin{1};
    temp = (mask).* double(preframeMov);
    area =  sum(sum(mask));
    trace = squeeze(double(sum(sum(temp))/area));
    
else
    trace = 0;
end
 



end