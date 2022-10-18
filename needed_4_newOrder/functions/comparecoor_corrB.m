function logicCorr = comparecoor_corrB(x1, y1, tolerance)
% modified from comparecoor_corr, GM, July 2020
keep = [];
trash = [];
%x1 = [Data.synapsis.positionCorrectedX];
%y1 = [Data.synapsis.positionCorrectedY];
x2 = x1;
y2 = y1;
logicCorr = true(1,length(x1));
k = 1;
dist = [];
for i = 1:length(x1)
    for j = 1:length(x2)
        transit(1,1)= x1(i);
        transit(1,2)= y1(i);
        transit(2,1)= x2(j);
        transit(2,2)= y2(j);
        if pdist(transit)<= 0 + tolerance && i ~= j
            dist(k,1) = i;
            dist(k,2) = j;
            %dist(k,3) = pdist(transit);
            k = k +1;
        end      
        
    end
end

for m = 1:length(dist)
    if  ismember( dist(m,1), trash) == false
        if ismember( dist(m,2), trash) == false
            trash = [trash dist(m,2)];
        end
        
        if ismember( dist(m,1), keep) == false
            keep = [keep dist(m,1)];
        end        
    end
end
logicCorr(trash) = false;
logicCorr = logicCorr';
%Data.synapsis = Data.synapsis(logicCorr);
% %% NAME AND SAVE
% if flagSave == 1 
%     fname = sprintf('spH_%s_%s', Data.InfoExp.date, Data.InfoExp.name);% filename
%     %fname = [fname '_2nd_C'];
%     name2C = 'Data';
%     %name2C = [Data.InfoExp.name '_2C'];
%     %eval([name2C,'= Data;']);% change the name of the variable
%     eval ([sprintf('save %s/%s.mat %s ''-v7.3''', strrep(Data.InfoExp.date, '_', '.'),fname,name2C)]);
% end
% 
% 


end