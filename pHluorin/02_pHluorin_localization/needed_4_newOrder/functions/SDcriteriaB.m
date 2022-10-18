%% Successful events were those that had amplitudes greater than twice the SD of X frames (⬃1.5 s) before the event.
%% SD of the whole traceDt


function locEv = SDcriteriaB (locEv, roi, flagDepurate)
x = 10;%10
%% remove evs before the first stim
locEv = locEv([locEv.frame] > 4);
synID = unique([locEv.synID]', 'stable');
dfALL = zeros(size(locEv))';
roiIndx = (1:1:length(roi))';

for syn = 1:length(synID)
    mskSyn = ismember([locEv.synID], synID(syn))';
    roimsk = ismember([roi.synID],synID(syn));
    roiIndxSyn = roiIndx(roimsk);
    trace = roi(roiIndxSyn).traceDt;    
    frames = [locEv(mskSyn).frame]';
    %SDtrace =  std(trace);    
    fSyn = trace(frames);
    SDtrace = zeros(length(frames),1);
    
    for i= 1:length(SDtrace)
        if frames(i) > (x + 1) %11
            SDtrace(i) = mean(trace((frames(i) - x):(frames(i) - 1))) + 2* std(trace((frames(i) - x):(frames(i) - 1)));
        else
            %maxVal = frames(i) - 1;
            SDtrace(i) = std(trace(1:(frames(i) - 1)));
        end
    end
    df = fSyn -SDtrace;
    dfALL(mskSyn) = df;
end
mskPositdf = dfALL > 0;

if flagDepurate == 1
    locEv = locEv(mskPositdf);
    %dfALL = dfALL(mskPositdf);
else
    for i = 1 : length(mskPositdf)
        locEv(i).SDb = mskPositdf(i);
    end
        
end

end 