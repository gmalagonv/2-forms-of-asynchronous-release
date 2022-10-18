function eqtrace = eqpreframes(trace,DurationPart,numAP)

%DurationGroupsBckg = 5; %% needed when localBckgFlag == 1 
numParts = numAP/(DurationPart-1); %% this is true if 1 Hz stim

% if DurationPart ~= 1
%     numParts = numAP/(DurationPart-1); %% this is true if 1 Hz stim
% else
%     
% end

%mean_noresp_dt = zeros(41,41,numBckg);
eqtrace = trace;
%fr = 1;
%indxLocalBckg =  ceil(frame/(DurationGroupsBckg*20));
%X = (1:1:length(trace));
for part = 1:numParts
    
    factor = 4 + (part - 1) * (DurationPart -1)* 20;
    %pntsIndx = (factor:20:(factor + (DurationPart - 1) * 20));
    axisX = [factor;(factor + (DurationPart - 1) * 20)];
    axisY = eqtrace(axisX);
    slope = (axisY(2) - axisY(1)) ./ (axisX(2) - axisX(1));
    b = axisY(1)-slope*axisX(1);
    
    pntsX = (factor+1:1:((factor + (DurationPart - 1) * 20)-1));
    pntsY =  eqtrace(pntsX)';
    %localBckgMov(:,:,(fr:localBckgIndx(end)+19)) = median(movie(:, :, localBckgIndx),3);
    %fr = localBckgIndx(end) + 20;
    lineSlope = slope.* pntsX + b;
    lineRef = ones(1,length(pntsY))* axisY(1);
    
    dif = pntsY - lineSlope;
    
    eqtrace(axisX) = axisY(1);
    eqtrace(pntsX) = lineRef + dif;
    eqtrace((axisX(2) + 1):end) = eqtrace((axisX(2) + 1):end) + (axisY(1) - axisY(2));
%     if gr ~= numParts
%         for i = fr: axisX(end) + 19
%             eqtrace(:,:,i) = median(movie(:, :, axisX),3);
%         end
%     else
%          for i = fr: length(movie)
%              eqtrace(:,:,i) = median(movie(:, :, axisX),3);
%          end
%     end
%     fr = axisX(end) + 20;

end

end






