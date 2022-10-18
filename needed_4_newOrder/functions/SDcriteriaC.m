%% Successful events were those that had amplitudes greater than twice the SD of the whole trace
%% SD from the 'pure noise' trace, meaning a trace composed by frames right before the stim  

function locEv = SDcriteriaC (locEv, roi, flagDepurate, InfoExp)
%% remove evs befroe the first stim
[~, firstPreframe] = max(InfoExp.descartFrame);
lastFrameBeforeTail = ((InfoExp.NumberAP)*((1/InfoExp.Frequency)/InfoExp.timeframe)) + firstPreframe;

AllFrame = (1:1:lastFrameBeforeTail)';
pureNoiseFrames = (AllFrame-(floor(AllFrame/20))*20)- firstPreframe;
pureNoiseFrames = AllFrame(pureNoiseFrames < 1);

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
    SDtrace =  std(trace(pureNoiseFrames));
    

    
%     trace = trace - min(trace) + 1; %% only positive values
%     for j = 1:length(preframesSynStart)
%         f0Syn(j,1) = mean(trace(preframesSynStart(j):preframesSynEnd(j)));    
%     end
    fSyn = trace(framesSyn);
    df = fSyn-(2*SDtrace);
    dfALL(mskSyn) = df;

end
mskPositdf = dfALL > 0;

if flagDepurate == 1
    locEv = locEv(mskPositdf);
    %dfALL = dfALL(mskPositdf);
else
    for i = 1 : length(mskPositdf)
        locEv(i).SDc = mskPositdf(i);
    end
        
end


end 