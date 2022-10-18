function all = corrDuplisynIDS(all)
if length([all.perSyn.synID])== length(unique([all.perSyn.synID])) && length([all.perSyn.synID]) == length(unique([all.allEv.synID]))
    disp ('no duplicated IDs')
    return
else

    %% correct duplicated synIDS:
    [~, indxUni] = unique([all.perSyn.synID], 'stable');
    indxSynsAll = 1:length(all.perSyn);
    dd_indx = setdiff(indxSynsAll, indxUni); %% dd_indx synapses with ID thta was already counted.
    ddID = unique([all.perSyn(dd_indx).synID], 'stable');

    %dd_indx = (dd_indx(dd1_indx'))';


    for dd = 1 : length(ddID)
        ID = ddID(dd);
        mskSynID = ismember([all.perSyn.synID], ID);
        indxSynsdd = indxSynsAll(mskSynID);
        %[~, ~, indexesON] = intervalsOnOff(ismember([all.perSyn.synID], ID));
        sameSyn = false;
      %  continSyn = false;
        if diff([all.perSyn(indxSynsdd).areaSynch]) == 0% meaning, if synapses come one after the other and have the same area. If so, probably both are the same synapse.

            sameSyn = true;
            all.perSyn(indxSynsdd) = [];

            if isfield(all, 'roi')
                all.roi(indxSynsdd) = [];
            end

            indxSynsAll = 1:length(all.perSyn);
        else
            for syndd = 1:length(indxSynsdd)
                newID = str2double([num2str(ID) num2str(syndd)]);
                all.perSyn(indxSynsdd(syndd)).synID = newID;
                if isfield(all, 'roi')
                    all.roi(indxSynsdd(syndd)).synID = newID;
                end
            end
        end



    %% correct synID allEv
        [~, ~, indexesON] = intervalsOnOff(ismember([all.allEv.synID], ID));


        for part = 1:length(indexesON)
            indxS = (indexesON{part})'; %% update INdx
            if sameSyn == true
                newID = 666666;
            else 
                newID = str2double([num2str(ID) num2str(part)]);
            end

            for p = 1:length(indxS)
                all.allEv(indxS(p)).synID = newID;
            end
        end

    end
    all.allEv([all.allEv.synID] == 666666) = [];


    if length([all.perSyn.synID])== length(unique([all.perSyn.synID])) && length([all.perSyn.synID]) == length(unique([all.allEv.synID]))
        disp('All duplicated synIDs corrected!')
    else

        disp('remaining duplicated synIDs')
        return
    end
end
end
