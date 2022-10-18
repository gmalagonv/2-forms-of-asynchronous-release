function DataorALL = distFromLastEv(DataorALL)

xy = [[DataorALL.allEv.x]',[DataorALL.allEv.y]'];
xyS= cat(1, [1 1], xy);
xyS(end, :) = [];

distFromLast = zeros(length(xyS) , 1);

for i = 1: length(distFromLast)
    distFromLast(i) = pdist(cat(1, xy(i,:), xyS(i,:)));
    
end


for syn = 1 :length(DataorALL.perSyn)
    [~, ~, indexesON]  = intervalsOnOff(ismember([DataorALL.allEv.synID], DataorALL.perSyn(syn).synID));
    for part = 1:length(indexesON)
        indxS = indexesON{part};
        distFromLast(indxS(1)) = NaN;
    end
end 

distFromLast = distFromLast * DataorALL.InfoExp.PixelSize;

for ev = 1:length(distFromLast)
    DataorALL.allEv(ev).distFromLast = distFromLast(ev);
end


end