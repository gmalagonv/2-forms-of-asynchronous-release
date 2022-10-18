%% Information
% load Synchronouse file "select folder"
% find localization of synapse from original 'ROI' and 'vGlut1' files
% making "ROI_Asynchro" file

%% Select folders 
close all;
clear;
set(0,'defaultAxesFontName', 'arial');

sel_fd_path = uigetdir('F:\DATA\Matlab\analysis');
cd(sel_fd_path);
thresh = inputdlg('What threshold do you want?', 'Threshold', [1, 35], {'10'});
thresh = thresh{1};
Name_ROI = dir('*ROI*.mat');
Name_ROI = {Name_ROI.name};
Name_ROI = Name_ROI(~contains(Name_ROI, 'Asynchro'));                                   % reject Asnchro files
Name_ROI = Name_ROI(~contains(Name_ROI, 'spo'));                                        % reject Asnchro files
Name_ROI = Name_ROI(~contains(Name_ROI, 'copy'));                                       % reject Asnchro files
Name_ROI = Name_ROI(contains(Name_ROI, ['thresh', thresh]));                            % reject Asnchro files

%% define variables
threshEv = 3.5;  
% decay values
decaytau = questdlg('What Fluoresence did you use?', 'Decay Tau', 'vGlut1 (330 ms)', 'syn-iGluSnFR (116 ms)', 'vGlut1 (330 ms)');
if strcmp(decaytau, 'syn-iGluSnFR (116 ms)')
    decaytau = 0.1162; 
elseif strcmp(decaytau, 'vGlut1 (330 ms)')
    decaytau = 0.330; 
end

%%
for i = 1:length(Name_ROI)
    % load ROI
    temp_ROI_name = strrep(Name_ROI{i}, '.', '_');                                      % 'ROI_01.31.2021_Data11_thresh10.mat' -> 'ROI_01_31_2021_Data11_thresh10'
    temp_ROI_name = erase(temp_ROI_name, '_mat');                                       % remove '.mat'
    load(Name_ROI{i});
    temp_ROI = eval(temp_ROI_name);                                                     % assign temp_ROI file
    eval(['clear ', temp_ROI_name]);
    WindowFrame = length(temp_ROI.InfoExp.descartFrame);
    % load vGlut1
    Name_vGlut1 = strrep(Name_ROI{i}, 'ROI', 'vGlut1-pHGFP');
    Name_vGlut1 = erase(Name_vGlut1, "_thresh" +  thresh);
    temp_vGlut1_name = strrep(Name_vGlut1, '.', '_');
    temp_vGlut1_name = strrep(temp_vGlut1_name, '-', '_');
    temp_vGlut1_name = erase(temp_vGlut1_name, '_mat');
    load(Name_vGlut1);
    temp_vGlut1 = eval(temp_vGlut1_name);                                               % assing temp_vGlut1 file 
    eval(['clear ', temp_vGlut1_name]);                                                 % erase Raw vGlut1_pHGFP_xx_xx_xxxx

    % checking decartframes
    if ~all(temp_ROI.InfoExp.descartFrame == [0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        if length(temp_ROI.InfoExp.descartFrame) ~= 20
            disp('Untached descart frames');
            return
        else
            disp(['descartFrame is ', num2str(temp_ROI.InfoExp.descartFrame)])
        end
    end
    
    % select individual synapse
    for ii = 1:temp_ROI.counted_sys_number
        LocTemp = cell(200, 1);
        Movie = double(temp_ROI.all(ii).movPerSyn);
        sigma = temp_ROI.InfoExp.sigma;
        index = temp_ROI.all(ii).indexData;
        
        % check whether correlation ROI vs vGlut1 file
        if ~all([temp_ROI.all(ii).positionCorrectedX, temp_ROI.all(ii).positionCorrectedY] ...
                == [temp_vGlut1.synapsis(index).positionCorrectedX, temp_vGlut1.synapsis(index).positionCorrectedY])
            return
        end

        MaskTemp = temp_vGlut1.synapsis(index).Mask;                        % mask
        % detection event
        for  iii = 1:temp_ROI.InfoExp.NumberAP
            % load movie and pre and post Frame
            postFrame = Movie(:, :, (iii - 1)*20 + 6);
            postFrame = postFrame + (1-exp(-50/decaytau))*(Movie(:, :, (iii - 1)*20 + 5) - mean(Movie(:, :, (iii - 1)*20 + 3:(iii - 1)*20 + 4), 3));        % compensation decay
            preFrame = Movie(:, :, (iii - 1)*20 + 5);
            [pstruct, ~, ~, ~] = FusionLoc(postFrame, preFrame, sigma, 'mask', double(MaskTemp), 'DeltaAcTresh', threshEv,'OnlyOne', true ); %-------------> DETECTION EVENT!!
            if ~isempty(pstruct)
                pstruct.EventN = iii;
                LocTemp{iii} = pstruct;
            end
        end
        % assign asynchronous information
        if isempty(LocTemp)                                                        % whether it has event
            Data.all(ii).positionCorrectedX = temp_ROI.all(ii).positionCorrectedX;
            Data.all(ii).positionCorrectedY = temp_ROI.all(ii).positionCorrectedY;
            Data.all(ii).NlocEvents = length(Loc.Event);                    % # of event
            Data.all(ii).events_localization = Loc;                         %
            Data.all(ii).indexData = temp_ROI.all(ii).indexData;
            Data.all(ii).roi_trace = temp_ROI.all(ii).roi_trace;
            Data.all(ii).roi_traceF = temp_ROI.all(ii).roi_traceF;
        else
            LocTemp = LocTemp(~cellfun('isempty',LocTemp));                         % rearrange LocTemp
            Loc.Event = cell2mat(LocTemp);                                          % Merge LocTemp
            Data.all(ii).positionCorrectedX = temp_ROI.all(ii).positionCorrectedX;
            Data.all(ii).positionCorrectedY = temp_ROI.all(ii).positionCorrectedY;
            Data.all(ii).NlocEvents = length(Loc.Event);                    % # of event
            Data.all(ii).events_localization = Loc;                         %
            Data.all(ii).indexData = temp_ROI.all(ii).indexData;
            Data.all(ii).roi_trace = temp_ROI.all(ii).roi_trace;
            Data.all(ii).roi_traceF = temp_ROI.all(ii).roi_traceF;
        end
        clear Loc
        % making trace (roi_traceF)
        for iii = 1:Data.all(ii).NlocEvents
            EventN = Data.all(ii).events_localization.Event(iii).EventN;
            Data.all(ii).trace(iii).FrameN = EventN;
            Data.all(ii).trace(iii).roi_traceF = Data.all(ii).roi_traceF((EventN-1) * WindowFrame + 1 : 1 : EventN * WindowFrame + 4);
            Data.all(ii).trace(iii).N_roi_traceF = Data.all(ii).trace(iii).roi_traceF/mean(Data.all(ii).trace(iii).roi_traceF([1:4]));
            Data.all(ii).trace(iii).roi_trace = Data.all(ii).roi_trace((EventN-1) * WindowFrame + 1 : 1 : EventN * WindowFrame + 4);
            Data.all(ii).trace(iii).N_roi_trace= Data.all(ii).trace(iii).roi_trace/mean(Data.all(ii).trace(iii).roi_trace([1:4]));
        end
    end
    Data.InfoExp = temp_ROI.InfoExp;
    Data = JY_map_Syn_Asyn(Data, Name_ROI{i});
    % save
    temp_ROI_name = strrep(temp_ROI_name, 'ROI', 'ROI_Asynchro');                                       % name 'ROI_Asynchro_01_31_2021_Data1_thresh10'
    eval([temp_ROI_name, ' = Data;']);
    Name_ROI_save = strrep(Name_ROI{i}, 'ROI', 'ROI_Asynchro');
    eval(['save(''', Name_ROI_save, ''', ''', temp_ROI_name, ''');']);
    clearvars -except Name_ROI sel_fd_path threshEv thresh decaytau
end
beep
pause(1)
beep
pause(1)
beep
