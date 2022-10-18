function [dfAll, msk, accFr] =  testSimpleDf(roi, lastFrameBeforeTail,allFlag)
%firstPreframe
x = 4;
accFr = cell(length(roi),1);
dfAll = zeros(length(roi),1);


if allFlag == 1
    %framesPost = ((firstPreframe + 1):1:lastFrameBeforeTail);
    framesPost = ((x + 1):1:lastFrameBeforeTail);
else
    %framesPost = ((firstPreframe + 1):20:lastFrameBeforeTail);
    framesPost = ((x + 1):20:lastFrameBeforeTail);
end
%


framesPreLast = framesPost -1;
%framesPreFirst = framesPreLast -(firstPreframe -1);
framesPreFirst = framesPreLast -(x -1);

for syn = 1:length(roi)
    trace = roi(syn).traceDt;
    f = trace(framesPost);
    f0 = zeros(size(f));
    for fr = 1: length(f)
        f0(fr) = mean(trace(framesPreFirst(fr):framesPreLast(fr)));
    end
        
    %f0 = trace(framesPreFirst:framesPreLast);
    df = (f - f0)./ f0;
%     if sum(df > 0) >= 70
%         msk(syn) = true;
%     else
%         msk(syn) = false;
%     end
    accFr{syn} =  framesPost(df>0);
    
    dfAll(syn) = sum(df>0)/length(framesPost);
end
msk = dfAll >= 0.5;

end