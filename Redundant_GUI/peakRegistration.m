function out = peakRegistration(synapticTraces, synSpace, rollWindow, threshFactor, peakInfluence)
    warning('off')
    
    % variables shared across nested functions
    contents.detectedPeaks = {}; 
    contents.registeredPeaks = {};
    currentSynapse = 1;
    frameLine = [];
    regPeaksIdxs = [];
    magentaTrace = [];
    
    % make figure
    screenSize = get(0, 'ScreenSize');
    regFig = uifigure('Units', 'pixels', 'Visible', 'on', ...
                'Position', [screenSize(3)/5 screenSize(4)/6 1300 790],... % [left bottom width height]
                'Name', 'Peak registration interface', 'ToolBar', 'none', 'MenuBar', 'none', 'NumberTitle', 'off');
    movegui(regFig, 'center'); 
    
    % plot all synapses in same figure
    axesAll = uiaxes(regFig, 'Position', [100 280 1100 500]);
    hold(axesAll, 'on');
     for ii = 1:size(synapticTraces,1)
         plot(axesAll, ones(size(synapticTraces(ii,:)))+(ii/synSpace)+synapticTraces(ii,:), 'k', 'linewidth', 0.8);
     end
     axesAll.YLabel.String = 'All synapses';
     set(axesAll, 'xtick', [], 'ytick', [], 'ylim', [0.5+(1/synSpace) 1.5+(ii/synSpace)],'xlim',[1 length(synapticTraces(1,:))])
    
    % plot first synapse
    axesOne = uiaxes(regFig, 'Position', [100 150 1100 120]);   
    updateTrace(currentSynapse)
    placeFrameLine(axesOne)
    
    % create interactive buttons
    frameSlider = uislider(regFig,'Position', [450 120 400 30],'limits', size(synapticTraces(1,:)), 'ValueChangingFcn', @(frameSlider, event) frameSliderMoved(frameSlider, event));
    prev = uibutton(regFig, 'push', 'position', [200 80 80 30], 'text', 'Prev', 'ButtonPushedFcn', @(prev, event) prevButtonPush(prev, event));
    next = uibutton(regFig, 'push', 'position', [300 80 80 30], 'text', 'Next', 'ButtonPushedFcn', @(next, event) nextButtonPush(next, event));
    peak = uibutton(regFig, 'push', 'position', [950 80 80 30], 'text', 'Peak', 'ButtonPushedFcn', @(peak, event) registerPeak(peak, event));
    submit = uibutton(regFig, 'push', 'position', [1050 80 80 30], 'text', 'Submit', 'ButtonPushedFcn', @(submit, event) submitPeaks(submit, event));
    textTitle = uilabel(regFig, 'Position',[200 40 200 20],'Text','Registered peaks for synapse: 1');
    textPeaks = uilabel(regFig, 'Position',[450 40 400 20],'Text','...');
    done = uibutton(regFig, 'push', 'position', [1050 40 80 30], 'text', 'Done', 'ButtonPushedFcn', @(done, event) finishAnalysis(done, event));

    % return output to workspace only when Done is pressed
    waitfor(done);
    if size(contents.registeredPeaks,2) ~= size(contents.detectedPeaks,2)
        addCells = abs(diff([size(contents.registeredPeaks,2),size(contents.detectedPeaks,2)]));
        contents.registeredPeaks{1,end+addCells} = zeros(size(synapticTraces(1,:)))';
    end
    for l = 1:size(contents.registeredPeaks,2)
        if isempty(contents.registeredPeaks{1,l}) 
            contents.registeredPeaks{1,l} = zeros(size(synapticTraces(1,:)))';
        end
    end
    out = contents;
    
    % nested functions
    function frameSliderMoved(frameSlider, event)
        if isempty(frameLine)
            placeFrameLine(axesOne)
        end
        frameLine.XData = [round(event.Value) round(event.Value)];
    end

    function prevButtonPush(prev, event)
        if currentSynapse ~= 1
            currentSynapse = currentSynapse - 1;
        end
        
        resetParams()
        updateTrace(currentSynapse)
        
        placeFrameLine(axesOne)
        textTitle.Text = sprintf('Registered peaks for synapse: %d', currentSynapse);
    end

    function nextButtonPush(next, event)
        if currentSynapse < size(synapticTraces,1)
            currentSynapse = currentSynapse + 1;
        end
        resetParams()
        updateTrace(currentSynapse)
        
        placeFrameLine(axesOne)
        textTitle.Text = sprintf('Registered peaks for synapse: %d', currentSynapse);
    end

    function updateTrace(synNum)
        plot(axesOne, synapticTraces(synNum, :),'k', 'linewidth', 1)
        hold(axesOne, 'on');
        contents.detectedPeaks{synNum} = getRasterActivity(synapticTraces(synNum, :), rollWindow, threshFactor, peakInfluence, 'values');
        [p, ~, v] = find(contents.detectedPeaks{synNum});
        plot(axesOne,p,v,'or')
        axesOne.XLabel.String = 'Frames';
        axesOne.YLabel.String = 'dF/F0';
        set(axesOne, 'ylim',[-0.1 0.5], 'xlim',[1 length(synapticTraces(1,:))])
        
        magentaTrace = plot(axesAll, ones(size(synapticTraces(synNum,:)))+(synNum/synSpace)+synapticTraces(synNum,:),'m', 'linewidth', 1);
    end

    function placeFrameLine(axObj)
        frameLine = plot(axObj, ones(2,1), axObj.YLim, 'b');
        frameSlider.Value = 1;
    end

    function registerPeak(peak, event)
        regPeaksIdxs(end+1) = frameLine.XData(1);
        textPeaks.Text = num2str(regPeaksIdxs);
    end

    function submitPeaks(submit, event)
        vals = synapticTraces(currentSynapse, regPeaksIdxs);
        pks = zeros(size(synapticTraces(currentSynapse,:)));
        pks(regPeaksIdxs) = vals;
        contents.registeredPeaks{currentSynapse} = pks';
    end

    function finishAnalysis(done, event)
        close(regFig)
    end

    function resetParams()
        hold(axesOne, 'off');
        magentaTrace.Color = 'k';
        magentaTrace.LineWidth = 0.8;
        textPeaks.Text = '';
        regPeaksIdxs = [];
    end
end

