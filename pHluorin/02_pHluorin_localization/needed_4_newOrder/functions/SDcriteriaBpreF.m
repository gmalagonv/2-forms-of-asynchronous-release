%% Successful events were those that had amplitudes greater than twice the SD of X frames (⬃1.5 s) before the event.
%% SD of the whole traceDt


function locEv = SDcriteriaBpreF (locEv, roi, flagDepurate)
x = 5;%10
%% remove evs before the first stim
locEv = locEv([locEv.frame] > 4);
synID = unique([locEv.synID]', 'stable');
logicMaskALL = true(size(locEv))';
roiIndx = (1:1:length(roi))';

for syn = 1:length(synID)
    mskSyn = ismember([locEv.synID], synID(syn))';
    roimsk = ismember([roi.synID],synID(syn));
    roiIndxSyn = roiIndx(roimsk);
    trace = roi(roiIndxSyn).traceDt;    
    frames = [locEv(mskSyn).frame]';
    frames = frames -1; % preframe
    %SDtrace =  std(trace);    
    fSyn = trace(frames);
    logicMask = true(length(frames),1);
    %SDtraceB = zeros(length(frames),1);
    
    for i= 1:length(logicMask)
        if frames(i) > (x + 1) %11
            %fSyn(i) = ((fSyn(i)-mean(trace((frames(i) - x):(frames(i) - 1))))^2)^(1/2);
            %SDtrace(i) = 2 * std(trace((frames(i) - x):(frames(i) - 1)));
            meanPiece = mean(trace((frames(i) - x):(frames(i) - 1)));
            stdPiece = std(trace((frames(i) - x):(frames(i) - 1)));
            
            %cond1 = mean(trace((frames(i) - x):(frames(i) - 1))) - std(trace((frames(i) - x):(frames(i) - 1)));
            
            %SDtrace(i) = mean(trace((frames(i) - x):(frames(i) - 1))) - std(trace((frames(i) - x):(frames(i) - 1)));
        else
           %fSyn(i)=((fSyn(i)-mean(trace(1:(frames(i) - 1))))^2)^(1/2);
            %SDtrace(i) = 2 * std(trace(1:(frames(i) - 1)));
            meanPiece = mean(trace(1:(frames(i) - 1)));
            stdPiece =  std(trace(1:(frames(i) - 1)));
            %SDtrace(i) = mean(trace(1:(frames(i) - 1))) - std(trace(1:(frames(i) - 1)));
        end
        cond1 = meanPiece - stdPiece;
        cond2 = meanPiece + stdPiece;
        logicMask(i) = fSyn(i) > cond1 && fSyn(i) < cond2; 
    end
    %df = fSyn-SDtrace;
    logicMaskALL(mskSyn) = logicMask;
end
%mskPositdf = logicMaskALL > 0;

if flagDepurate == 1
    locEv = locEv(logicMaskALL);
    %dfALL = dfALL(mskPositdf);
else
    for j = 1 : length(logicMaskALL)
        locEv(j).SDpreF = logicMaskALL(j);
    end
        
end

end 