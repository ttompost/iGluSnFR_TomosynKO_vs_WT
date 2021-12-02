classdef SniffApp < matlab.apps.AppBase
    properties (Access = public)
        Figure
        FileMenu
        FileMenuOpen
        AnalysisMenu
        RegisterAnalysis
        SaveAnalysis
        AxesFov
        AxesMov
        AxesRoi
        FilesOpen
        ListFiles
        RoiSlider
        FrameSlider
        MinRoiSlider
        MaxRoiSlider
        MinFrameSlider
        MaxFrameSlider
        AxesTraces
        SynapseNum
        ShowingSynapse
        RemoveRoi
        AddMean
        ShowAllRois
        ShowAllRoiTraces
        AnalysisMethodsButtons
        GradientMethod
        GradientMethodApply
        ThresholdMethod
        ThresholdMethodApply
        FovTitle
        RecTitle
        FluoTitle
        RoiTitle
    end
    
    properties (Access = private)
        filesDir   % directory with all files  (fov, rois, mov, traces)
        filesDataStore
        dataTable
        fovFile     % name of the FOV file
        fovFileOld  % to check if new file will be loaded
        roiFile     % directory path to the rois folder
        movFile     % name of the recording file
        movFileOld
        movStack    % recording
        allRois     % names and fiji-generated coordinates
        nRois
        nFrames
        fovRaw          % raw FOV
        fovEdited       % raw FOV edited
        fovOverlayed    % fovEdited overlayed with Fiji-detected rois
        roiTraces       % gray pixel values
        roiFzeros       % f0 for each roi
        roiDeltaFsRaw   % dF/f0 for each roi
        timeLine        % moving line to follow with frame slider
        selectedRoi     % roi shown on the FOV/still image
        timeVec         % time vector
        newPosition     % new roi selection using slider
        selectedFile    % file to be shown in gui
        recFrequency
        recExposure
        threshResults
        roiAnalysisResults
        roiTraceOld
    end
    
    % callback methods
    methods (Access = private)
        function FileMenuOpenSelected(sapp, event)
            %%% this function asks the user to provide directory containing
            %%% files and displays all recordings (.tif), still images
            %%% (.tif and/or .jpg), fiji-traces (.csv) and fiji-rois
            %%% (RoiSet folder name)
            sapp.filesDir = uigetdir(pwd, 'Select data folder');
            addpath(genpath(sapp.filesDir));
            cd(sapp.filesDir)
            set(0, 'DefaulttextInterpreter', 'none');
            sapp.filesDataStore = imageDatastore (sapp.filesDir, 'FileExtensions', {'.tif','.jpg','.csv'},'IncludeSubfolders', false);
            [~, fileNames, fileExt] = cellfun(@fileparts, sapp.filesDataStore.Files, 'UniformOutput', false);
            
            % some file names start with "._", so here we delete them
            if any(cell2mat(cellfun(@(x) startsWith(x, '._'), fileNames, 'uniformoutput', false))) == 1
                [namesWithPrefixes, ~] = find(cell2mat(cellfun(@(x) startsWith(x, '._'), fileNames, 'uniformoutput', false)) == 1); 
                for nameWPref_idx = namesWithPrefixes
                    fileNames(nameWPref_idx) = [];
                    fileExt(nameWPref_idx) = [];
                end
            end
            fullfileNames = cellfun(@(x, y) join([char(x), char(y)]), fileNames, fileExt, 'uniformoutput', false);
            
            % imageDatastore does not list folder names, so here we search for a folder with rois
            if ismac || isunix 
                roisFolder = dir('**/RoiSet*'); 
            elseif ispc
                roisFolder = dir('**\RoiSet*');
            end
            
            if ~isempty(roisFolder)
                for roiSubs = 1:size(roisFolder,1)
                    fullfileNames{end+1} = roisFolder(roiSubs).name;
                end
            end
            sapp.ListFiles.String = fullfileNames;
            sapp.FilesOpen.Enable = 'on';
        end
        
        function OpenFilesPress(sapp, event)
            %%% opens any selected name from the list and displays it in
            %%% the predetermined axes
            if sapp.FilesOpen.Value == 1
                sapp.selectedFile = char(sapp.ListFiles.String(sapp.ListFiles.Value));
                fileTypeMsg = 'This is a:'; %sprintf('What type of file is \n%s? \n', sapp.selectedFile);
                typeOptions = {'FOV/still image', 'Recording', 'Roi set', 'Roi traces/Results'};
                fileType = listdlg('PromptString', fileTypeMsg, 'SelectionMode','single','ListString',typeOptions);
%                 idx = 1;
%                 if startsWith(sapp.selectedFile, 'Results_') || startsWith(sapp.selectedFile, 'RoiSet_')
%                     idx=idx+1;
%                 end
%                 nameComponents = split(sapp.selectedFile, '_');
%                 recWeek = weeknum(datetime(nameComponents{idx},'inputformat','yyMMdd'));
%                 recID = nameComponents{idx+1};
%                 csID = nameComponents{idx+2};
%                 fovID = split(nameComponents{idx+3},'.');
%                 fovID = fovID{1};
%                 sapp.dataTable = table(sapp.selectedFile, typeOptions(fileType), recWeek, csID, recID, fovID);
%                 sapp.dataTable.Properties.VariableNames = {'FileName', 'FileType', 'WeekID', 'CoverslipID', 'RecID', 'FovID'};
                if fileType == 1 % fov
                    if isempty(sapp.fovFileOld)
                        openNewFov = 1;
                    elseif strcmp(sapp.selectedFile, sapp.fovFileOld) == 0
                        dMessage = 'Pressing OK will change the current FOV and clear rest of the axes. Continue?';
                        dTitle = 'Confirm change of FOV';
                        selection = questdlg(dMessage, dTitle,'OK','Cancel','Cancel'); % last one is a default option
                        if strcmp(selection,'OK')  == 1
                            openNewFov = 1;
                            axes(sapp.AxesFov); cla;
                            axes(sapp.AxesMov); cla;
                            axes(sapp.AxesTraces); cla;
                            axes(sapp.AxesRoi); cla;
                            sapp.RoiSlider.Enable = 'off';
                            sapp.FrameSlider.Enable = 'off';
                            sapp.MaxRoiSlider.String = 'x';
                            set(sapp.ShowAllRois, 'Enable', 'off', 'Value', 0);
                            set(sapp.AddMean, 'value',0, 'enable', 'off');
                            set(sapp.ShowAllRoiTraces, 'value', 0, 'enable', 'off');
                            sapp.fovOverlayed = [];
                            sapp.allRois = [];
                            sapp.movStack = [];
                            sapp.roiTraces = [];
                            sapp.roiDeltaFsRaw = [];
                            sapp.timeLine = [];
                        else
                            openNewFov = 0;
                        end
                    else
                        openNewFov = 0;
                    end
                    if openNewFov == 1
                        loadTime = waitbar(0, sprintf('Loading: %s', sapp.selectedFile));
                        sapp.fovFileOld = sapp.selectedFile;
                        sapp.fovRaw = imread(sapp.selectedFile);
                        if length(size(sapp.fovRaw)) == 2 % grayscale (or single channel)
                            sapp.fovRaw = im2uint8(sapp.fovRaw);
                            fovTophat = imtophat(sapp.fovRaw,strel('square',7)); 
                            sapp.fovEdited = imadjust(imcomplement(fovTophat) - sapp.fovRaw);
                        else
                            sapp.fovEdited = sapp.fovRaw;
                        end
                        imshow(sapp.fovEdited, 'Parent', sapp.AxesFov);
                        delete(loadTime);
                        set(sapp.SynapseNum, 'String', 'none');
                    else
                        axes(sapp.AxesFov);
                        imshow(sapp.fovEdited);
                        if ~isempty(sapp.allRois)
                            thisRoi = get(sapp.RoiSlider, 'value');
                            sapp.selectedRoi = rectangle('position',[sapp.allRois(thisRoi).coords(2), sapp.allRois(thisRoi).coords(1),...
                                abs(sapp.allRois(thisRoi).coords(1)-sapp.allRois(thisRoi).coords(3)), abs(sapp.allRois(thisRoi).coords(2)-sapp.allRois(thisRoi).coords(4))],...
                                'edgecolor','r','linewidth',1.2);
                        else
                            set(sapp.RoiSlider,'Enable','off');
                            set(sapp.SynapseNum, 'String', 'none');
                        end
                    end

                % open results
                elseif fileType == 4
                    if isempty(sapp.roiTraces) || strcmp(sapp.selectedFile, sapp.roiTraceOld) == 0
                        sapp.roiTraceOld = sapp.selectedFile;
                        axes(sapp.AxesTraces); cla;
                        if ismac || isunix
                            sapp.roiTraces = readtable(join([sapp.filesDir '/' sapp.selectedFile]));
                        elseif ispc
                            sapp.roiTraces = readtable(join([sapp.filesDir '\' sapp.selectedFile]));
                        end
                        loadTime = parfor_wait(sapp.nFrames, 'Waitbar', true,'FileName','screen');
                        sapp.roiTraces = sapp.roiTraces{:,2:end}'; % rois are rows now, frames are columns
                        sapp.roiTraces = flip(sapp.roiTraces); % apparently, fiji stores the list in reverse?
                        if isempty(sapp.movStack) % if traces are loaded before recording
                            sapp.nFrames = size(sapp.roiTraces,2); 
                            sapp.nRois = size(sapp.roiTraces,1);
                        end
                        if sapp.nRois ~= size(sapp.roiTraces,1) && sapp.nFrames ~= size(sapp.roiTraces,2)
                            return
                        else % calculate df/f0 from raw traces
                            seg=10; % 10 segments
                            segBreaks=[1:round(sapp.nFrames/seg):sapp.nFrames, sapp.nFrames];
                            sapp.roiFzeros=zeros(sapp.nRois,seg); % nRoi x segments
                            for segNum=1:seg
                                sapp.roiFzeros(:,segNum)=min(sapp.roiTraces(:,segBreaks(segNum):segBreaks(segNum+1)-1),[],2);
                            end
                            sapp.roiFzeros=repmat(mean(sapp.roiFzeros,2),1,sapp.nFrames);
                            sapp.roiDeltaFsRaw=zeros(size(sapp.roiTraces));     % dF/F0 for all synapses
                            axes(sapp.AxesTraces); 

                            for roiNum=1:sapp.nRois
                                dF=(sapp.roiTraces(roiNum,:)-sapp.roiFzeros(roiNum,:))./ sapp.roiFzeros(roiNum,:);
                                sapp.roiDeltaFsRaw(roiNum,:)=dF;

                                plot(dF,'color',[0.63 0.63 0.63]);  hold on;
                                loadTime.Send;
                            end
                            plot(mean(sapp.roiDeltaFsRaw,1),'color','r');
                        end
                        loadTime.Destroy;
                    else
                        axes(sapp.AxesTraces); 

                        for roiNum=1:sapp.nRois
                            plot(sapp.roiDeltaFsRaw(roiNum,:),'color',[0.63 0.63 0.63]); hold on;
                        end
                        if get(sapp.AddMean, 'value') == 1
                            plot(mean(sapp.roiDeltaFsRaw,1),'color','r');
                        end
                        axes(sapp.AxesRoi); cla;
                    end
                    if ~isempty(sapp.movStack)
                        sapp.timeLine.YData = sapp.AxesTraces.YLim;
                        uistack(sapp.timeLine, 'top');
                    end
                    set(sapp.SynapseNum,'string', 'all');
                    editTracesAxesLabels(sapp);

                    axes(sapp.AxesRoi); cla;
                    set(sapp.AddMean, 'value',1, 'enable', 'on');
                    set(sapp.ShowAllRoiTraces, 'value', 1, 'enable', 'on');
                    placeFrameTracer(sapp);
                        
                    % open recording    
                elseif fileType == 2
                    sapp.movFile = sapp.selectedFile;
                    recParamsQuestions = {'Imaging frequency (in Hz):','Exposure time (in ms):'};
                    recParamsTitle = 'Recording parameters';
                    recParamsDefInput = {'10', '0.100'}; % 10Hz, 100ms exposure 
                    recParamsUserInput = inputdlg(recParamsQuestions, recParamsTitle, [1 36], recParamsDefInput);
                    sapp.recFrequency = str2double(recParamsUserInput{1});
                    sapp.recExposure = str2double(recParamsUserInput{2});
                    if isempty(sapp.movStack) || strcmp(sapp.movFile, sapp.movFileOld) == 0
                        sapp.movFileOld = sapp.movFile;
                        movInfo = imfinfo(sapp.movFile);
                        sapp.nFrames = length(movInfo);
                        loadTime = parfor_wait(sapp.nFrames, 'Waitbar', true,'FileName','screen');
                        sapp.timeVec = 0:sapp.recExposure:(sapp.nFrames/sapp.recFrequency); 
                        recTime = seconds(sapp.nFrames/sapp.recFrequency);
                        recTime.Format = 'mm:ss';
                        parfor frame = 1 : sapp.nFrames
                            mov(:,:,:,frame) = imresize(imread(sapp.movFile, frame),[size(sapp.fovRaw,1), size(sapp.fovRaw,2)]);
                            loadTime.Send;
                        end
                        pools = gcp('nocreate'); % get currently running parallel pools
                        delete(pools); % shut the pools down
                        sapp.movStack = mov; clear mov;
                        imshow(sapp.movStack(:,:,:,1),'Parent',sapp.AxesMov);
                        loadTime.Destroy;

                        set(sapp.FrameSlider,'Enable','on','Min',1,'Max', sapp.nFrames,'Value',1,...
                        'SliderStep', [1/(sapp.nFrames-1), 10/(sapp.nFrames-1)]);
                        set(sapp.MaxFrameSlider,'String',join([string(recTime) 'min']));
                    else
                        imshow(sapp.movStack(:,:,:,1),'Parent',sapp.AxesMov);
                        sapp.FrameSlider.Value = 1;
                    end

                    if isempty(sapp.timeLine)
                        placeFrameTracer(sapp);
                    end
                
                % rois folder
                elseif fileType == 3 % if selection points to a roiset folder
                    if isempty(sapp.fovOverlayed)
                        if ismac || isunix
                            sapp.roiFile = dir(join([sapp.filesDir '/' sapp.selectedFile]));
                        elseif ispc
                            sapp.roiFile = dir(join([sapp.filesDir '\' sapp.selectedFile]));
                        end
                        loadTime = waitbar(0, sprintf('Loading: %s', sapp.selectedFile));
                        sapp.fovOverlayed = sapp.fovEdited;
                        sapp.allRois = struct();
                        sapp.nRois = size(sapp.roiFile,1)-2;
                        for roiNum = 3:size(sapp.roiFile,1)
                            if startsWith(sapp.roiFile(roiNum).name, '._')
                                continue
                            else
                                roi = ReadImageJROI(sapp.roiFile(roiNum).name);
                                sapp.allRois(roiNum-2).name = roi.strName; % roi name generated in fiji
                                sapp.allRois(roiNum-2).coords = roi.vnRectBounds; % fiji coordinates

                                sapp.fovOverlayed = insertShape(sapp.fovOverlayed,'rectangle',[roi.vnRectBounds(2), roi.vnRectBounds(1),...
                                    abs(roi.vnRectBounds(1)-roi.vnRectBounds(3)), abs(roi.vnRectBounds(2)-roi.vnRectBounds(4))],...
                                    'color','blue');
                                waitbar(roiNum/sapp.nRois,loadTime);
                            end
                        end
                        set(sapp.RoiSlider,'Enable','on','Min',1,'Max',sapp.nRois,'Value',1,...
                            'SliderStep', [1/(sapp.nRois-1), 10/(sapp.nRois-1)]);
                        set(sapp.MaxRoiSlider,'String',num2str(sapp.nRois));
                        delete(loadTime);
                    end
                    imshow(sapp.fovOverlayed,'Parent',sapp.AxesFov);

                    axes(sapp.AxesFov);
                    highlightRoi(sapp);
                    sapp.SynapseNum.String = sapp.RoiSlider.Value;
                    set(sapp.ShowAllRois, 'Enable', 'on', 'Value', 1);
                end
            end
        end
        
        function RegisterAnalysisSelected(sapp, event)
            %%%
            thisRoi = sapp.RoiSlider.Value;
            thisFile = split(sapp.fovFileOld,'.');
            gatheredResults = {thisFile{1}, sapp.allRois(thisRoi).name, sapp.allRois(thisRoi).coords,sapp.roiDeltaFsRaw(thisRoi,:)};
            generalData = cell2table(gatheredResults, 'VariableNames', {'FileName','RoiID', 'RoiLocation', 'dFoverF'});
            if isempty(sapp.roiAnalysisResults) % create new
                sapp.roiAnalysisResults = [generalData sapp.threshResults];
            else % append
                sapp.roiAnalysisResults = [sapp.roiAnalysisResults; [generalData sapp.threshResults]];
            end
        end
        
        function SaveAnalysisResultsSelected(sapp, event)
            %%% 
            fileName = split(sapp.fovFileOld, '.');
            fileName = join(['GuiResults_' fileName{1}]);
            [file,path] = uiputfile(fileName);
        end
        
        function RoiSliderMoved(sapp, event)
            %%% slider movement runs through individual synapses and plots
            %%% its df/f trace
            if sapp.RoiSlider.Visible
                if nargin == 2 && isprop(event, 'HorizontalScrollCount')
                    thisRoi = round(sapp.RoiSlider.Value + event.HorizontalScrollCount);
                    if thisRoi < 1
                        return
                    end
                else
                    thisRoi = round(sapp.RoiSlider.Value);
                end
                sapp.RoiSlider.Value = thisRoi;
                if ~isempty(sapp.roiDeltaFsRaw)
                    %set(sapp.RemoveRoi, 'enable', 'on');
                    axes(sapp.AxesTraces); hold off;
                    plot(sapp.roiDeltaFsRaw(thisRoi,:),'k','linewidth',1.2);
                    hold on;
                    if get(sapp.AddMean, 'value') == 1
                        plot(mean(sapp.roiDeltaFsRaw,1),'color','r');
                    end
                    editTracesAxesLabels(sapp);
                    placeFrameTracer(sapp);
                end
                sapp.newPosition = [sapp.allRois(thisRoi).coords(2), sapp.allRois(thisRoi).coords(1),...
                    abs(sapp.allRois(thisRoi).coords(1)-sapp.allRois(thisRoi).coords(3)), abs(sapp.allRois(thisRoi).coords(2)-sapp.allRois(thisRoi).coords(4))];
                sapp.selectedRoi.Position = sapp.newPosition;
                sapp.selectedRoi.EdgeColor = 'r';
                sapp.SynapseNum.String = thisRoi;
                sapp.ShowAllRoiTraces.Value = 0;
                sapp.GradientMethodApply.Enable = 'on';
                sapp.ThresholdMethodApply.Enable = 'on';
                
                % show zoomed in roi
                if ~isempty(sapp.movStack)
                    zoomRoi = imcrop(sapp.movStack(:,:,:,1),sapp.newPosition);
                    imshow(zoomRoi,'Parent',sapp.AxesRoi);
                end
            end
        end
        
        function FrameSliderMoved(sapp, event)
             %%% slider movement displays different recording frames 
             if sapp.FrameSlider.Visible 
                 if nargin == 2 && isprop(event, 'HorizontalScrollCount')
                    frame = round(sapp.FrameSlider.Value + event.HorizontalScrollCount);
                    if frame < 1
                        return
                    end
                else
                    frame = round(sapp.FrameSlider.Value);
                 end
                sapp.FrameSlider.Value = frame;
                if ~isempty(sapp.movStack)
                    axes(sapp.AxesMov);
                    imshow(sapp.movStack(:,:,:,frame));%,'Parent',sapp.AxesMov);
                    if length(sapp.AxesRoi.Children) == 1
                        axes(sapp.AxesRoi);
                        imshow(imcrop(sapp.movStack(:,:,:,frame),sapp.newPosition));%,'Parent',sapp.AxesRoi);
                    end
                end
                sapp.timeLine.XData = [frame frame];
             end
        end
        
        function GradientPress(sapp, event)
            %%% detects peaks by filtering the trace using gradient
            %%% approach and findpeaks function
            if sapp.GradientMethodApply.Value == 1
                % ask for user input
                gradParamsQuestions = {'Imaging frequency (in Hz):','Minimal prominence factor (MPF in: (Mean+(Std*MPF)):'};
                gradParamsTitle = 'Gradient method parameters';
                gradParamsDefInput = {num2str(sapp.recFrequency), '1.5'}; 
                gradParamsUserInput = inputdlg(gradParamsQuestions, gradParamsTitle, [1 54], gradParamsDefInput);
                
                % work with user input
                fps = str2double(gradParamsUserInput{1});
                minPromFac = str2double(gradParamsUserInput{2});
                synapse = str2double(get(sapp.SynapseNum,'String'));
                trace = sapp.roiDeltaFsRaw(synapse,:);
                grad = gradient(gradient(trace));
                filtered = (trace.*grad).*-1;
                axes(sapp.AxesTraces); hold off;
                [pks, idxs] = findpeaks(filtered,fps,'minpeakheight', mean(filtered)+(std(filtered)*minPromFac));
                plot(filtered,'-k', 'linewidth', 1.2); hold on; 
                %plot(trace,'color',ones(3,1).*0.63,'linewidth',2); % y axis will get too big
                if get(sapp.AddMean, 'value') == 1
                        plot(mean(sapp.roiDeltaFsRaw,1),'color','r');
                end
                plot(idxs*fps, pks, '*r','markersize',7);
                editTracesAxesLabels(sapp);
                placeFrameTracer(sapp);
            end
        end
        
        function ThresholdPress(sapp, event)
            %%% finds peaks based on the dynamic threshold method
            if sapp.ThresholdMethodApply.Value == 1
                % ask for user input
                threshParamsQuestions = {'Smoothing window lag (in ms):','Peak detection threshold (a.u.):', 'New peak influence (a.u.):'};
                threshParamsTitle = 'Threshold method parameters';
                treshParamsDefInput = {'50', '3.5', '0.2'}; 
                treshParamsUserInput = inputdlg(threshParamsQuestions, threshParamsTitle, [1 54], treshParamsDefInput);
                
                % work with user input
                windowLag = str2double(treshParamsUserInput{1});
                thresh = str2double(treshParamsUserInput{2});
                influence = str2double(treshParamsUserInput{3});
                synapse = str2double(get(sapp.SynapseNum,'String'));
                trace = sapp.roiDeltaFsRaw(synapse,:);
                [pks, avgFilt, stdFilt] = ThresholdingAlgo(trace, windowLag, thresh, influence);
                
                % gather params
                tResults = {windowLag, thresh, influence, pks};
                sapp.threshResults = cell2table(tResults, 'VariableNames',{'threshParam_Lag','threshParam_Threshold',...
                    'threshParam_Influence','threshParam_Peaks'});
                
                axes(sapp.AxesTraces); hold off;
                plot(trace,'-k','linewidth', 1.2); hold on;
                if get(sapp.AddMean, 'value') == 1
                        plot(mean(sapp.roiDeltaFsRaw,1),'color','r');
                end
                for posPk = find(pks==1)
                    plot(posPk,trace(posPk),'*r','markersize',7); 
                end
                ylim(sapp.AxesTraces.YLim);
                editTracesAxesLabels(sapp);
                placeFrameTracer(sapp);
            end
        end
        
%         function RemoveRoiPress(sapp, event)
%             %%%
%         end
        
        function AddMeanPress(sapp, event)
            %%% adds or removes mean from all roi traces
            axes(sapp.AxesTraces)
            switch sapp.AddMean.Value
                case 1
                    plot(mean(sapp.roiDeltaFsRaw,1),'color','r');
                case 0
                   plotContent = sapp.AxesTraces.Children;
                   for plt = 1:length(plotContent)
                       if isequal(plotContent(plt).Color, [1 0 0]) && strcmp(plotContent(plt).LineStyle, '-') % searching for a red line (i.e., mean)
                           delete(plotContent(plt))
                       end
                   end
            end
        end
        
        function AllRoisPress(sapp, event)
            %%% adds or removes all rois from fov axes
            axes(sapp.AxesFov)
            switch sapp.ShowAllRois.Value
                case 1
                    imshow(sapp.fovOverlayed,'Parent',sapp.AxesFov);
                    highlightRoi(sapp);
                case 0
                    imshow(sapp.fovEdited,'Parent',sapp.AxesFov);
                    highlightRoi(sapp);
            end
        end
        
        function AllTracesPress(sapp, event)
            %%%
            axes(sapp.AxesTraces)
            switch sapp.ShowAllRoiTraces.Value
                case 1
                    hold off
                    for roiNum=1:sapp.nRois
                        hold on
                        plot(sapp.roiDeltaFsRaw(roiNum,:),'color',[0.63 0.63 0.63]);
                    end       
                    if get(sapp.AddMean, 'value') == 1
                        plot(mean(sapp.roiDeltaFsRaw,1),'color','r');
                    end
                    ylim([min(min(sapp.roiDeltaFsRaw)) max(max(sapp.roiDeltaFsRaw))]);
                    editTracesAxesLabels(sapp);
                    placeFrameTracer(sapp);
                case 0
                    axes(sapp.AxesTraces); hold off;
                    plot(sapp.roiDeltaFsRaw(sapp.RoiSlider.Value,:),'color','k','linewidth',1.2);
                    sapp.SynapseNum.String = num2str(sapp.RoiSlider.Value);
                    hold on;
                    if get(sapp.AddMean, 'value') == 1
                        plot(mean(sapp.roiDeltaFsRaw,1),'color','r');
                    end
                    editTracesAxesLabels(sapp);
                    placeFrameTracer(sapp);
                    sapp.GradientMethodApply.Enable = 'on';
                    sapp.ThresholdMethodApply.Enable = 'on';
            end
        end
        function highlightRoi(sapp)
            %%% highlights currently selected roi on fov in red
            axes(sapp.AxesFov);
            thisRoi = sapp.RoiSlider.Value;
            sapp.selectedRoi = rectangle('position',[sapp.allRois(thisRoi).coords(2), sapp.allRois(thisRoi).coords(1),...
                        abs(sapp.allRois(thisRoi).coords(1)-sapp.allRois(thisRoi).coords(3)), abs(sapp.allRois(thisRoi).coords(2)-sapp.allRois(thisRoi).coords(4))],...
                        'edgecolor','r','linewidth',1.2);
        end
        
        function placeFrameTracer(sapp)
            %%% places a vertical bar to follow the frame slider
            ylims=sapp.AxesTraces.YLim;
            sapp.timeLine = plot(sapp.AxesTraces, ones(2,1)*sapp.FrameSlider.Value, ylims, 'b');
            set(sapp.AxesTraces,'ylim',ylims);
        end
        
        function editTracesAxesLabels(sapp)
            %%% takes care of labeling the xy axes for roi traces
            sapp.AxesTraces.YLabel.String = '\DeltaF/F_0';
            sapp.AxesTraces.YLabel.FontSize = 14;
            sapp.AxesTraces.YLabel.Interpreter = 'tex';
            sapp.AxesTraces.XLabel.String = 'Time (s)';
            sapp.AxesTraces.XLabel.FontSize = 14;
%            set(gca,'XTickLabel',sapp.timeVec(1:sapp.nFrames/sapp.recFrequency:end),'fontsize',14);
        end
      
        function makeGui(sapp)
            %%% creates all components in SniffApp
            
            % define figure size
            screenSize = get(0, 'ScreenSize');
            sapp.Figure = figure('Units', 'pixels', 'Visible', 'on', ...
                'Position', [screenSize(3)/5 screenSize(4)/6 1300 790],... % [left bottom width height]
                'Name', 'iGluSnFR data visualisation app', 'ToolBar', 'none', 'MenuBar', 'none', 'NumberTitle', 'off');
            set(sapp.Figure, 'WindowStyle','normal'); movegui(sapp.Figure, 'center');  
            % create the menu bar
            sapp.FileMenu = uimenu(sapp.Figure, 'Text', 'File');
            sapp.FileMenuOpen = uimenu(sapp.FileMenu, 'Text', 'Open data folder',...
                'MenuSelectedFcn', createCallbackFcn(sapp, @FileMenuOpenSelected, true));
            sapp.AnalysisMenu = uimenu(sapp.Figure, 'Text', 'Analysis');
            sapp.RegisterAnalysis = uimenu(sapp.AnalysisMenu, 'Text', 'Register results',...
                'MenuSelectedFcn', createCallbackFcn(sapp, @RegisterAnalysisSelected, true));
            sapp.SaveAnalysis = uimenu(sapp.AnalysisMenu, 'Text', 'Save results',...
                'MenuSelectedFcn', createCallbackFcn(sapp, @SaveAnalysisResultsSelected, true));
            % place image axes 
            sapp.AxesFov = axes('Units', 'pixels', 'Position', [30 340 400 400],'Visible','on'); 
            sapp.AxesFov.XTick = [];
            sapp.AxesFov.YTick = [];
            sapp.AxesMov = axes('Units', 'pixels', 'Position', [455 340 400 400],'Visible','on');
            sapp.AxesMov.XTick = [];
            sapp.AxesMov.YTick = [];
            sapp.AxesRoi = axes('Units', 'pixels', 'Position', [905 470 70 70],'Visible','on');
            sapp.AxesRoi.XTick = [];
            sapp.AxesRoi.YTick = [];
            sapp.AxesTraces = axes('Units', 'pixels', 'Position', [80 70 1180 150],'Visible','on');
            % place files boxes
            sapp.ListFiles = uicontrol(sapp.Figure, 'Style', 'listbox','Units', 'pixels', 'Position', [890 600 360 170]);
            % place interactive/open buttons
            sapp.FilesOpen = uicontrol(sapp.Figure, 'Style', 'pushbutton',...
                'String', 'Open', 'Units', 'pixels', 'Position',[1180 560 70 35],...
                'FontSize', 12, 'Enable', 'off','Callback', createCallbackFcn(sapp, @OpenFilesPress, true));
            % place sliders
            sapp.FrameSlider = uicontrol(sapp.Figure, 'Style', 'slider','Units', 'pixels',...
                'Position', [495 290 320 30],'Visible', 'on','Enable','off',...
                'Callback', createCallbackFcn(sapp, @FrameSliderMoved, true));
            addlistener(sapp.FrameSlider,'Value','PreSet',@(~,~)FrameSliderMoved(sapp));
            sapp.RoiSlider = uicontrol(sapp.Figure, 'Style', 'slider','Units', 'pixels',...
                'Position', [70 290 320 30],'Visible', 'on','Enable','off',...
                'Callback', createCallbackFcn(sapp, @RoiSliderMoved, true));
            addlistener(sapp.RoiSlider,'Value','PreSet',@(~,~)RoiSliderMoved(sapp));
            % place uneditable slider text boxes
            sapp.MinFrameSlider = uicontrol(sapp.Figure, 'Style', 'text','String', '0:00',...
                'Units', 'pixels', 'Position', [455 290 30 30],'FontSize', 12, 'Enable', 'on');
            sapp.MaxFrameSlider = uicontrol(sapp.Figure, 'Style', 'text','String', 'x',...
                'Units', 'pixels', 'Position', [820 290 50 30],'FontSize', 12, 'Enable', 'on');
            sapp.MinRoiSlider = uicontrol(sapp.Figure, 'Style', 'text','String', '1',...
                'Units', 'pixels', 'Position', [30 290 30 30],'FontSize', 12, 'Enable', 'on');
            sapp.MaxRoiSlider = uicontrol(sapp.Figure, 'Style', 'text','String', 'x',...
                'Units', 'pixels', 'Position', [400 290 30 30],'FontSize', 12, 'Enable', 'on');
            sapp.ShowingSynapse = uicontrol(sapp.Figure, 'Style', 'text','String', 'Showing synapse:',...
                'Units', 'pixels', 'Position', [160 270 100 30],'FontSize', 12, 'Enable', 'on');
            sapp.SynapseNum = uicontrol(sapp.Figure, 'Style', 'text','String', 'x',...
                'Units', 'pixels', 'Position', [270 270 30 30],'FontSize', 12, 'Enable', 'on');
            % titles
            sapp.FovTitle = uicontrol(sapp.Figure, 'Style', 'text','String', 'Field of view',...
                'Units', 'pixels', 'Position', [90 750 300 20],'FontSize', 14,'fontweight','bold', 'Enable', 'on');
            sapp.RecTitle = uicontrol(sapp.Figure, 'Style', 'text','String', 'Recording',...
                'Units', 'pixels', 'Position', [510 750 300 20],'FontSize', 14,'fontweight','bold', 'Enable', 'on');
            sapp.FluoTitle = uicontrol(sapp.Figure, 'Style', 'text','String', 'Fluorescence trace(s)',...
                'Units', 'pixels', 'Position', [520 230 300 30],'FontSize', 16,'fontweight','bold', 'Enable', 'on');
            sapp.RoiTitle = uicontrol(sapp.Figure, 'Style', 'text','String', 'Selected ROI',...
                'Units', 'pixels', 'Position', [870 555 150 20],'FontSize', 14,'fontweight','bold', 'Enable', 'on');
            % button controls
            sapp.RemoveRoi = uicontrol(sapp.Figure, 'Style', 'radiobutton',...
                'String', 'Exclude this ROI', 'Units', 'pixels', 'Position', [1010 510 200 30],...
                'FontSize', 12, 'Enable', 'off','callback', createCallbackFcn(sapp, @RemoveRoiPress, true));
            sapp.ShowAllRois = uicontrol(sapp.Figure, 'Style', 'radiobutton',...
                'String', 'Show all ROIs', 'Units', 'pixels', 'Position', [1010 480 200 30],...
                'FontSize', 12, 'Enable', 'off','callback',createCallbackFcn(sapp, @AllRoisPress, true));
            sapp.AddMean = uicontrol(sapp.Figure, 'Style', 'radiobutton',...
                'String', 'Show all-ROIs mean', 'Units', 'pixels', 'Position', [950 250 200 30],...
                'FontSize', 12, 'Enable', 'off','callback',createCallbackFcn(sapp, @AddMeanPress, true));
            sapp.ShowAllRoiTraces = uicontrol(sapp.Figure, 'Style', 'radiobutton',...
                'String', 'Show all ROI traces', 'Units', 'pixels', 'Position', [950 280 200 30],...
                'FontSize', 12, 'Enable', 'off','callback',createCallbackFcn(sapp, @AllTracesPress, true));
            % analysis buttons
            sapp.AnalysisMethodsButtons = uibuttongroup(sapp.Figure, 'Units', 'pixels',...
                'Title', 'Trace analysis', 'Position', [890 330 350 120], 'FontSize', 12);
            sapp.GradientMethod = uicontrol(sapp.Figure, 'Style', 'text','String','Gradient filtering',...
                'Units', 'pixels', 'Position', [920 390 200 30],'FontSize', 12, 'Enable', 'on','HorizontalAlignment','left');
            sapp.GradientMethodApply = uicontrol(sapp.Figure, 'Style', 'pushbutton',...
                 'String', 'Use', 'Units', 'pixels', 'Position', [1150 392 50 30],...
                 'FontSize', 12, 'Enable', 'off','Callback', createCallbackFcn(sapp, @GradientPress, true));
            sapp.ThresholdMethod = uicontrol(sapp.Figure, 'Style', 'text','String','Threshold detection',...
                'Units', 'pixels', 'Position', [920 350 200 30],'FontSize', 12, 'Enable', 'on','HorizontalAlignment','left');
            sapp.ThresholdMethodApply = uicontrol(sapp.Figure, 'Style', 'pushbutton',...
                 'String', 'Use', 'Units', 'pixels', 'Position', [1150 352 50 30],...
                 'FontSize', 12, 'Enable', 'off','Callback', createCallbackFcn(sapp, @ThresholdPress, true));
        end
    end
    
    methods (Access = public)
        function sapp = SniffApp
            %%% creates the SniffApp
            makeGui(sapp);
        end
    end
    
end
