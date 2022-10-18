function randVal = randFromdistribution(NumEvs, distEvs2sim)

[density, bins] = histcounts(distEvs2sim);
bins2 = (bins(1:(length(bins)-1)) + bins(2:length(bins)))./2;
val = randsample(bins2,NumEvs,true,density);

%% refine values
p = @(r) find(r<bins,1,'first'); % function
indxBins = arrayfun(p,val);

upLim = bins(indxBins);
lowLim = bins(indxBins-1);

randVal = zeros(length(upLim), 1);
% 
% a = 50;
% b = 100;
% r = (b-a).*rand(1000,1) + a;


for i = 1:length(randVal)
    randVal(i) = (upLim(i)- lowLim(i)) .* rand(1,1) + lowLim(i);
end


end