function DataorALL = timeFromLastEv(DataorALL)

frames = [DataorALL.allEv.frame]';
framesS = cat(1, 1, frames);
framesS(end) = [];
timeFromLast = frames- framesS;


for syn = 1 :length(DataorALL.perSyn)
    [~, ~, indexesON]  = intervalsOnOff(ismember([DataorALL.allEv.synID], DataorALL.perSyn(syn).synID));
    for part = 1:length(indexesON)
        indxS = indexesON{part};
        timeFromLast(indxS(1)) = NaN;
    end
end 

timeFromLast = timeFromLast * 50;

for ev = 1:length(timeFromLast)
    DataorALL.allEv(ev).timeFromLast = timeFromLast(ev);
end


end