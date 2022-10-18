function temp_ROI_Asyn = JY_map_Syn_Asyn(temp_ROI_Asyn, Name_Exp_file_mat)
%% finding Asynchrnouse events
% Using ROI_Asynchro to find synapse
% load StuffLoc
StuffLoc_file = strrep(Name_Exp_file_mat, 'ROI', 'StuffLoc');
load(StuffLoc_file); %#ok<LOAD>
temp_StuffLoc_File = erase(StuffLoc_file, '.mat');
temp_StuffLoc_File = strrep(temp_StuffLoc_File, '.', '_');
temp_StuffLoc = eval(temp_StuffLoc_File);
PixelSize = temp_ROI_Asyn.InfoExp.PixelSize;

for ii = 1:length(temp_StuffLoc.perSyn)
    index = temp_StuffLoc.perSyn(ii).synNum;
    Syn_Loc_XY = temp_StuffLoc.perSyn(ii).XY * PixelSize;                                                   % Synchronouse events localization
    Syn_cluster_XY = reshape([temp_StuffLoc.perSyn(ii).cluster.centroid], 2, [])' * PixelSize;
    centerPhyXY = [temp_StuffLoc.perSyn(ii).centerPhyX; temp_StuffLoc.perSyn(ii).centerPhyY]' * PixelSize;   % center of AZ
    area_point = temp_StuffLoc.perSyn(ii).area_point * PixelSize;
    index = [temp_ROI_Asyn.all.indexData] == index;                                             %finding identical synapse ex) [0, 0, 1, 0, 0,...]
    temp_ROI_Asyn.all(index).Syn_Loc_XY = Syn_Loc_XY;
    temp_ROI_Asyn.all(index).Syn_Cluster_XY = Syn_cluster_XY;
    temp_ROI_Asyn.all(index).centerPhyXY = centerPhyXY;
    temp_ROI_Asyn.all(index).area_point = area_point;
    Syn_dist = pdist([centerPhyXY; Syn_Loc_XY]);
    temp_ROI_Asyn.all(index).Syn_dist = Syn_dist(1:size(Syn_Loc_XY, 1));                                                                        % Distance center to Synchronouse 
    if nnz(index) == 1 && temp_ROI_Asyn.all(index).NlocEvents ~= 0                                                                             % Only exist Asynchronouse events
        Asyn_loc_XY = [[temp_ROI_Asyn.all(index).events_localization.Event.x]',...                  % Asynchronouse events localization
            [temp_ROI_Asyn.all(index).events_localization.Event.y]'] * PixelSize;
        for iii = 1:size(Asyn_loc_XY, 1)                                                            % check each asynchronouse event
            %             [crossx, crossy] = polyxpoly([area_point(:, 1); area_point(1, 1)], [area_point(:, 2); area_point(1, 2)], ...
            %                 [centerPhyXY(1, 1), Asyn_loc_XY(iii, 1)], [centerPhyXY(1, 2), Asyn_loc_XY(iii, 2)]);
            [Norm_Fac, ectopic] = func_Norm_Fac(Asyn_loc_XY(iii, :), centerPhyXY, area_point);
            %             Norm_Fac = pdist([centerPhyXY; crossx, crossy]) + 25;                       % Normalization fctor
            if  ectopic
                temp_ROI_Asyn.all(index).events_localization.Event(iii).ectopic = ectopic;
                temp_ROI_Asyn.all(index).events_localization.Event(iii).dist = pdist([centerPhyXY; Asyn_loc_XY(iii, :)]);
                temp_ROI_Asyn.all(index).events_localization.Event(iii).Norm_Fac = Norm_Fac;
                temp_ROI_Asyn.all(index).events_localization.Event(iii).Norm_dist = pdist([centerPhyXY; Asyn_loc_XY(iii, :)])/Norm_Fac;
            else
                temp_ROI_Asyn.all(index).events_localization.Event(iii).ectopic = ectopic;
                temp_ROI_Asyn.all(index).events_localization.Event(iii).dist = pdist([centerPhyXY; Asyn_loc_XY(iii, :)]);
                temp_ROI_Asyn.all(index).events_localization.Event(iii).Norm_Fac = Norm_Fac;
                temp_ROI_Asyn.all(index).events_localization.Event(iii).Norm_dist = pdist([centerPhyXY; Asyn_loc_XY(iii, :)])/Norm_Fac;
            end
        end
    end


end

