function Data = addFieldSyn2dEvs(Data, field)
indxAllEv =(1:length(Data.allEv)); 
for syn = 1:length(Data.perSyn)
    val = Data.perSyn(syn).(field);
    mskSyn = ismember([Data.allEv.synID], Data.perSyn(syn).synID);
    indxSynEv = indxAllEv(mskSyn);
    
    for ev = 1:length(indxSynEv)
        Data.allEv(indxSynEv(ev)).(field) = val;
        
    end
    

end
