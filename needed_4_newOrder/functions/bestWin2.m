function locEv = bestWin2 (locEv, minNumSynchEv, windowDuration, InfoExp, flagDepurate)

[~, firstFrameStim] = max(InfoExp.descartFrame);
firstFrameStim = firstFrameStim + 1;
NumframPerStim = (1/InfoExp.Frequency)/InfoExp.timeframe;

windowLenght = windowDuration * InfoExp.Frequency;
windowNum = InfoExp.NumberAP - windowLenght + 1; 


synID = unique([locEv.synID]', 'stable');
framesALL = [locEv.frame]';
synchALL = [locEv.synch]';




mskFinalFrames = false(length(framesALL),1); 
%mskFinalSyn = true(length(synID),1);



for syn = 1:length(synID)
    mskSyn = ismember([locEv.synID], synID(syn))';
    %maskLocalFrames = true(sum(mskSyn), 1);
    frameSyn = framesALL(mskSyn);
    
    %mskFrms = false(frameSyn);
    synchFrames = frameSyn(synchALL(mskSyn));
    synchStims = ((synchFrames - firstFrameStim)/NumframPerStim) + 1;
    onesZeros = zeros(InfoExp.NumberAP,1);
    onesZeros(synchStims) = 1;
    windowScore = zeros(windowNum, 1);
    
    for i = 1:windowNum
        windowScore(i,1) = sum(onesZeros(i: i + windowLenght-1));
    end
    
    [~, maxScoreLoc] = max(windowScore);
    limitStims = [maxScoreLoc, (maxScoreLoc + windowLenght -1)];
    limitFrames = ((limitStims - 1) * NumframPerStim) + firstFrameStim;
    limitFrames(2) = limitFrames(2) + NumframPerStim -1; %% All the frames after the event
    mskFrms = (frameSyn >= limitFrames(1) & frameSyn <= limitFrames(2));
    if sum(mskFrms(synchALL(mskSyn))) >= minNumSynchEv
        mskFinalFrames(mskSyn) = mskFrms;
    end
    
    
end
if flagDepurate == 1
    locEv = locEv(mskFinalFrames);
else
    for p = 1: length(mskFinalFrames)
        locEv(p).bestWin = mskFinalFrames(p);
    end
end


end 