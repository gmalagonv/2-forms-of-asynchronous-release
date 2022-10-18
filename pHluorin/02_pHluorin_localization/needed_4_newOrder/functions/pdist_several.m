function dist = pdist_several(x1,y1,x2, y2)
%function dist = comparecoor(Data1,Data2, tolerance)
%x1 = [Data1.synapsis.positionCorrectedX];
%y1 = [Data1.synapsis.positionCorrectedY];
%x2 = [Data2.synapsis.positionCorrectedX];
%y2 = [Data2.synapsis.positionCorrectedY];
k = 1;
dist = [];
for i = 1:length(x1)
    for j = 1:length(x2)
        transit(1,1)= x1(i);
        transit(1,2)= y1(i);
        transit(2,1)= x2(j);
        transit(2,2)= y2(j);        
        %if pdist(transit)<= 0 + tolerance %&& i ~= j
        dist(k,1) = i;
        dist(k,2) = j;
        dist(k,3) = pdist(transit);
        k = k +1;   
      
        
    end
end



end