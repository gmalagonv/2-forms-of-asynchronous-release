%% Data load from Asynchronouse ROI data

%% pre setup
close all; set(0, 'defaultAxesFontName', 'arial');

%% load summary file 
ComputerName = getComputerName;                                                                             % Computer name
answer_original = questdlg("Do you want to analize from Datadeloaded Analized data");
if strcmp(answer_original, 'No')

elseif strcmp(answer_original, 'Yes')                                                                       % in case 
    if strcmp(ComputerName, 'jongyunmyeongdesktop')
        cd('E:\Google Drive\making paper\Synapse\temporal_analysis');                                       % change CD
        [Name_Exp_file_mat, path] = uigetfile;
        cd(path); 
        load(Name_Exp_file_mat);
        Name_Exp_file = erase(Name_Exp_file_mat,  '.mat');
%         load(Name_Exp_file);
        date = eval([Name_Exp_file, '.info.date']);
        for i = 1:length(date)                                                                              % load ROI_Asynchro
            file_name = strrep(date{i}, 'ROI', 'ROI_Asynchro');
            temp_day = char(extractBetween(file_name, 'ROI_Asynchro_', '_Data'));                         % temp_day is xx_xx_xxxx
            day = strrep(temp_day, '_', '.');                                                             % day is xx.xx.xxxx
            day_split = split(day, '.'); year = day_split{3};
            cd(['F:\DATA\Matlab\analysis\', year, '\', day]);
            load(strrep(file_name, temp_day, day))
        end
    else
       disp('It is not Jongyun''s destktop') 
        return
    end
elseif strcmp(answer_origincal, 'Cancel')
    return
end
ectopic = [];
Norm_dist = [];
ectopic_dist = [];
inside_dist = [];
Syn_dist = [];
% Syn_Loc_XY = [];
% ASyn_Loc_XY = [];
ROI_Asyn_list = who('ROI_Asynchro_*');
fg = figure(3);
ax = axes(fg);
map_answer = questdlg('Would you like to make a map?', 'AZ Map', 'Yes', 'No', 'No');                        % Asking whether 
n = 1;
for i = 1:size(ROI_Asyn_list, 1)
    temp_Asyn = eval(ROI_Asyn_list{i});
    for ii = 1:length(temp_Asyn.all)
        if temp_Asyn.all(ii).NlocEvents ~= 0  && any(strcmp(fieldnames(temp_Asyn.all(ii).events_localization.Event), 'ectopic'))                                                              % Igonore zero asynchronouse event 
            ectopic = [ectopic, [temp_Asyn.all(ii).events_localization.Event.ectopic]];                     %#ok<AGROW> <------- gethering ectopic event
            Norm_dist = [Norm_dist, [temp_Asyn.all(ii).events_localization.Event.Norm_dist]];               %#ok<AGROW> <------- gethering Norm_dist
            ect_idx = [temp_Asyn.all(ii).events_localization.Event.ectopic];
            ectopic_dist = [ectopic_dist, [temp_Asyn.all(ii).events_localization.Event(ect_idx).dist]];     %#ok<AGROW> 
            inside_dist = [inside_dist, [temp_Asyn.all(ii).events_localization.Event(~ect_idx).dist]];      %#ok<AGROW> 
            Syn_dist = [Syn_dist, temp_Asyn.all(ii).Syn_dist];                                              %#ok<AGROW>
            if strcmp(map_answer, 'Yes') && size(temp_Asyn.all(ii).Syn_Loc_XY, 1) > 10 && temp_Asyn.all(ii).NlocEvents > 5
                pgon = polyshape(temp_Asyn.all(ii).area_point(:, 1), temp_Asyn.all(ii).area_point(:, 2));
                plot(pgon, 'FaceColor',  [255, 236, 191]/255, 'EdgeColor', [255, 212, 127]/255);
                ax.XAxis.Visible = 'off';
                ax.YAxis.Visible = 'off';
                set(ax, 'Color', 'none');
                hold(ax, 'on');
                scatter(ax, temp_Asyn.all(ii).Syn_Loc_XY(:, 1), temp_Asyn.all(ii).Syn_Loc_XY(:, 2), 40, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0, 0, 0]);
                scatter(ax, [temp_Asyn.all(ii).events_localization.Event.x]*temp_Asyn.InfoExp.PixelSize, ...
                    [temp_Asyn.all(ii).events_localization.Event.y]*temp_Asyn.InfoExp.PixelSize ...
                    , 40, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0, 0.6, 0]);
                scatter(ax, temp_Asyn.all(ii).Syn_Cluster_XY(:, 1), temp_Asyn.all(ii).Syn_Cluster_XY(:, 2), 40, 'Marker', 'o', 'MarkerEdgeColor', 'k');
                scatter(ax, temp_Asyn.all(ii).Syn_Cluster_XY(:, 1), temp_Asyn.all(ii).Syn_Cluster_XY(:, 2), 40, 'Marker', 'x', 'MarkerEdgeColor', 'k');
                scatter(ax, temp_Asyn.all(ii).centerPhyXY(1), temp_Asyn.all(ii).centerPhyXY(2), 40, 'Marker', 'd', 'markerEdgecolor', 'k', 'markerfacecolor', [255, 85, 0]/255);
                pause;
                hold(ax, 'off');
            end
            %             for iii = 1:temp_Asyn.all(iii).NlocEvents
%                 
%             end
            n = n + 1;
        end
    end
end

fg_1 = figure(1);
set(fg_1, 'units', 'inch', 'position', [1, 1,3.2, 2.2], 'color', 'w');
ax_1 = axes(fg_1);
Norm_dist_cdf = cdfplot(Norm_dist);                                 %, 'Color', [0, 153/255, 0])
Norm_dist_cdf.Color = [0, 0.6, 0];
Norm_dist_cdf.LineWidth = 2;
set(ax_1, 'unit', 'inch', 'GridColor', 'none', 'xlim', [0, 4], 'xtick', 0:4, 'ytick', 0:0.2:1, 'fontsize', 10, 'box', 'off', 'color', 'w');
ax_1.XLabel.String = 'Norm dist. to the AZ center';
ax_1.Title.Visible = 'off';
ax_1.YLabel.String = 'cumulative probability';
ax_1.YTickLabel{1} = '0.0';
ax_1.YTickLabel{end} = '1.0';
ax_1.TickDir = 'out';
ax_1.TickLength = [0.01, 0.5];
ax_1.Position([3, 4]) = [2.2857, 1.6];

fg_2 = figure(2);
set(fg_2, 'unit', 'inch', 'position', [1, 5, 3.2, 2.2], 'color', 'w');
ax_2 = axes(fg_2);
Syn_dist_cdf = cdfplot(Syn_dist);
Syn_dist_cdf.Color = [0, 0, 0];
Syn_dist_cdf.LineWidth = 2;
hold(ax_2, 'on');
Inside_dist_cdf = cdfplot(inside_dist);
Inside_dist_cdf.Color = [0, 0.6, 0];
Inside_dist_cdf.LineWidth = 2;
set(ax_2, 'unit', 'inch', 'GridColor', 'none', 'xlim', [0, 400], 'xtick', 0:100:400, 'ylim', [-0.01, 1], 'ytick', 0:0.2:1, 'fontsize', 10, 'box', 'off', 'color', 'none');
ax_2.XLabel.String = 'dist. to the AX center (nm)';
ax_2.Title.Visible = 'off';
ax_2.YLabel.String = 'cumulative probability';
ax_2.YTickLabel{1} = '0.0';
ax_2.YTickLabel{end}  = '1.0';
ax_2.TickDir = 'out';
ax_2.TickLength = [0.01, 0.5];
ax_2.Position([3, 4]) = [2.2857, 1.6];

% [a, b] = kstest2(Syn_dist, inside_dist);