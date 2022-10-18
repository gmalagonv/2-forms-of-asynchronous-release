%% Successful events were those that had amplitudes greater than twice the SD of the whole trace
%% SD of the whole traceDt

function locEv = SDcriteria (locEv, roi, flagDepurate)
%% remove evs befroe the first stim
locEv = locEv([locEv.frame] > 4);

synID = unique([locEv.synID]', 'stable');
dfALL = zeros(size(locEv))';
roiIndx = (1:1:length(roi))';

for syn = 1:length(synID)
    mskSyn = ismember([locEv.synID], synID(syn))';
    framesSyn = [locEv(mskSyn).frame]';
%     preframesSynStart = framesSyn -4;
%     preframesSynEnd = framesSyn -1;
%     f0Syn = zeros(length(preframesSynStart),1);
    roimsk = ismember([roi.synID],synID(syn));
    roiIndxSyn = roiIndx(roimsk);
    trace = roi(roiIndxSyn).traceDt;
    SDtrace = mean(trace(5:end)) + 2*std(trace(5:end));
    
%     trace = trace - min(trace) + 1; %% only positive values
%     for j = 1:length(preframesSynStart)
%         f0Syn(j,1) = mean(trace(preframesSynStart(j):preframesSynEnd(j)));    
%     end
    fSyn = trace(framesSyn);
    df = fSyn-SDtrace;
    dfALL(mskSyn) = df;

end
mskPositdf = dfALL > 0;

if flagDepurate == 1
    locEv = locEv(mskPositdf);
    %dfALL = dfALL(mskPositdf);
else
    for i = 1 : length(mskPositdf)
        locEv(i).SDa = mskPositdf(i);
    end
        
end


end 