function Data = addSomeFitVariables(Data)
%% add the density of events per synapse
for j = 1:length(Data.perSyn)
    Data.perSyn(j).synchEvDensity = Data.perSyn(j).NumSynch/Data.perSyn(j).areaSynch; 
    
end

% for i = 1:length(Data.allEv)
%     Data.allEv(i).valminusc = Data.allEv(i).val - Data.allEv(i).c;
%     if  Data.allEv(i).valminusc < 0
%         Data.allEv(i).valminusc = NaN;
%     end
%     Data.allEv(i).errFit = sqrt((Data.allEv(i).val - (Data.allEv(i).A + Data.allEv(i).c))^2);
%     Data.allEv(i).errFit2 = Data.allEv(i).errFit/ Data.allEv(i).val;
%     Data.allEv(i).dA = Data.allEv(i).A / Data.allEv(i).c;
% 
%     
%     
% end