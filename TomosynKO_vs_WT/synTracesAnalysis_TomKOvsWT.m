%% add paths
yourPath = '';
addpath(genpath(yourPath));

%% load main data
clear all; close all;

load TomosynKO_RawData.mat

tomVSwt = detrendSynapticTraces(tomVSwt, 'tomVSwt',[99, 101]); % windows were chosen based on a priori exploratory analysis
tomVSwt = getEvokedAmplitudes(tomVSwt, 'tomVSwt');

%% Plot mean traces per condition
f=figure('name','Mean fluorescence traces'); set(f,'windowstyle','docked');
subplot(1,4,[1 2])
stdShade(cell2mat(table2cell(tomVSwt(tomVSwt.Condition == '225','dFoverF0_detrend'))),0.1,[0 0 1]); hold on
stdShade(cell2mat(table2cell(tomVSwt(tomVSwt.Condition == '148','dFoverF0_detrend'))),0.1,[0 0 0])
legend({'tomosynKO std','tomosynKO mean','wild-type std','wild-type mean'},'location','northwest');
xlim([1 800])
xlabel('frames (1frame=0.1s)'); ylabel('iGluSnFR dF/F0');

subplot(1,4,3)
stdShade(cell2mat(table2cell(tomVSwt(tomVSwt.Condition == '225','dFoverF0_detrend'))),0.1,[0 0 1]); hold on
stdShade(cell2mat(table2cell(tomVSwt(tomVSwt.Condition == '148','dFoverF0_detrend'))),0.1,[0 0 0])
xlim([801 1010])
xa = get(gca,'xaxis');
xlabel('frames (1frame=0.1s)'); ylabel('iGluSnFR dF/F0');

subplot(1,4,4)
stdShade(cell2mat(table2cell(tomVSwt(tomVSwt.Condition == '225','dFoverF0_detrend'))),0.1,[0 0 1]); hold on
stdShade(cell2mat(table2cell(tomVSwt(tomVSwt.Condition == '148','dFoverF0_detrend'))),0.1,[0 0 0])
xlim([1011 1110])
xlabel('frames (1frame=0.1s)'); ylabel('iGluSnFR dF/F0');

%% Plot boxplots with weekly datapoints for 0.75 and 2 Hz stimulations
plotTypes = {'fiveAt075Hz_dFoF0','fiveAt2Hz_dFoF0'};
genotypes = [148, 225]; colors =  ['k';'b'];
weeks = [22 25 26];

for pt = 1:length(plotTypes)
%     wt=mean(mean(cell2mat(table2cell(tomVSwt(tomVSwt.Condition == '148',plotTypes(pt)))))); % mean value from all wt measurements,
%     potentially to be used for normalization
    f=figure('name',string(plotTypes(pt))); clf; set(f,'windowstyle','docked');
    for g = 1:length(genotypes)
        allMeans = [];
        %clf
        subplot(121+g-1); hold on;
        for w = 1:length(weeks)
            weekData = table2cell(tomVSwt(tomVSwt.Condition == num2str(genotypes(g)) & tomVSwt.Week == weeks(w),plotTypes(pt)));
            for ap=1:5
                means =  cellfun(@(x) mean(x(:,ap)), weekData,'uniformoutput', false);
                stds = cellfun(@(x) std(x(:,ap)), weekData,'uniformoutput', false);

                if w == 1; shift = -0.2; shape = 's'; elseif w == 3; shift = 0.2; shape = 'o'; else; shift = 0; shape =  'd'; end

                plot(ones([1 length(means)])+shift+ap-1,cell2mat(means),[shape colors(g)],'markersize',7);
                %errorbar(ones([1 length(means)])+shift+ap-1,cell2mat(means),cell2mat(stds),color);

                if ap == 1
                    lastSyn = size(allMeans,1);
                end
                allMeans(lastSyn+1:lastSyn+length(means),ap) = cell2mat(means);
            end
            %return
        end
        if genotypes(g) == 148; title('wild-type'); elseif genotypes(g) == 225; title('tomosynKO');end
        xlabel('AP #'); ylabel('iGluSnFR dF/F0');
        boxplot(allMeans,[1:1:5],'notch','on','PlotStyle','traditional','colors',colors(g),'Widths',0.9,'MedianStyle','line');
        %ylim([-0.05 1]) % same y axes for better comparison, but then some outliers are cropped
        set(gca,'yminortick','on');
    end
end

%% Calculate synaptic participation per cell
% load the already existing data:
load TomosynKO_SynapticParticipation_075Hz.mat % participation during 0.75 Hz
load TomosynKO_SynapticParticipation_2Hz.mat % participation during 2 Hz
    % and add this data to the table:
tomVSwt = [tomVSwt, tomVswt_raster075data, tomVswt_raster2data];

%%%  this code runs the calculation from scratch (about 30min with 2 workers)
% stimStarts_075 = [600 615 630 645 660];
% frontTail_075 = 4;
% backTail_075 = 9;
% minScore_075 = 7;
% 
% stimStart_2 = 810;
% frontTail_2 = 3;
% backTail_2 = 4;
% minScore_2 = 5;

% parfor cellNum = 1:size(tomVSwt,1)
%     apScores_075 = []; apScores_2 = [];
%     apPositions_075 = []; apPositions_2 = [];
%     apBinary_075 = []; apBinary_2 = [];
%     for synNum = 1:size(tomVSwt.dFoverF0_detrend{cellNum,1},1)
%         for apTime = 1:length(stimStarts_075)
%             noiseStart075 = randi(580);
%             noiseTrace075 = tomVSwt.dFoverF0_detrend{cellNum,1}(synNum,noiseStart075:noiseStart075+backTail_075);
%             [apBin_075, apScore_075, apPos_075] = isItAP(tomVSwt.dFoverF0_detrend{cellNum,1}(synNum,stimStarts_075(apTime)-frontTail_075:stimStarts_075(apTime)+backTail_075),...
%                 noiseTrace075, 'ap', minScore_075); 
%             
%             apScores_075(synNum, apTime) = apScore_075;
%             if apPos_075 ~= 0
%                apPos_075 = apPos_075 + (stimStarts_075(apTime)-frontTail_075-1);
%             end 
%             apPositions_075(synNum, apTime) = apPos_075;
%             apBinary_075(synNum, apTime) = apBin_075;
%             
%             % analyse only the first AP for 2Hz
%             if apTime == 1 
%                 noiseStart2 = randi([880 980]);
%                 noiseTrace2 = tomVSwt.dFoverF0_detrend{cellNum,1}(synNum,noiseStart2:noiseStart2+backTail_2+frontTail_2);
%                 [apBin_2, apScore_2, apPos_2] = isItAP(tomVSwt.dFoverF0_detrend{cellNum,1}(synNum,stimStart_2-frontTail_2:stimStart_2+backTail_2),...
%                     noiseTrace2, 'ap', minScore_2); 
%                 apScores_2(synNum, apTime) = apScore_2;
%                 if apPos_2 ~= 0
%                     apPos_2 = apPos_2 + (stimStart_2-frontTail_2-1);
%                 end 
%                 apPositions_2(synNum, apTime) = apPos_2;
%                 apBinary_2(synNum, apTime) = apBin_2;
%             end
%         end
%     end
%     allApEvals_bin_075{cellNum,1} = apBinary_075;
%     allApEvals_scores_075{cellNum,1} = apScores_075;
%     allApEvals_pos_075{cellNum,1} = apPositions_075;
%     
%     allApEvals_bin_2{cellNum,1} = apBinary_2;
%     allApEvals_scores_2{cellNum,1} = apScores_2;
%     allApEvals_pos_2{cellNum,1} = apPositions_2;
% end
% try % if ran from scratch
%     tomVSwt = [tomVSwt, table(allApEvals_bin_075, allApEvals_scores_075, allApEvals_pos_075, allApEvals_bin_2, allApEvals_scores_2, allApEvals_pos_2, ...
%         'variablenames',{'binaryRaster075','apScores075','apPositions075','binaryRaster2','apScores2','apPositions2'})];
% catch ME
% end
% fprintf('Done.\n')

%% Plot synaptic participation per cell
genotypes = [148, 225];
colors =  ['k';'b'];
f=figure('name','Synaptic participation'); clf; hold on; set(f,'windowstyle','docked'); 
for g = 1:length(genotypes)
    fracts075 = cell2mat(cellfun(@(x) sum(x,1)./size(x,1), table2cell(tomVSwt(tomVSwt.Condition == num2str(genotypes(g)),'binaryRaster075')), 'uniformoutput', false));
    fracts2 = cell2mat(cellfun(@(x) sum(x,1)./size(x,1), table2cell(tomVSwt(tomVSwt.Condition == num2str(genotypes(g)),'binaryRaster2')), 'uniformoutput', false));
    s=subplot(131+g-1); hold on;
    for ap = 1:5
        subplot(s);
        shift = ap;
        jit = randi([-300 300]+(shift*1000), [1,length(fracts075(:,ap))])/1000;
        for syn = 1:length(fracts075(:,ap))
            scatter(jit(syn), fracts075(syn,ap),36,colors(g))
        end
        plot([-0.4 0.4]+shift, [median(fracts075(:,ap)) median(fracts075(:,ap))],['-' colors(g)]); hold on;
        
        if ap==1 % plot only 1st AP  for 2Hz
            subplot(133); hold on;
            shift = ap+g;
            jit = randi([-300 300]+(shift*1000), [1,length(fracts2(:,ap))])/1000;
            for syn = 1:length(fracts2(:,ap))
                scatter(jit(syn), fracts2(syn,ap),36,colors(g))
            end
            plot([-0.4 0.4]+shift, [median(fracts2(:,ap)) median(fracts2(:,ap))],['-' colors(g)]); hold on;
            title('2Hz')
            ya = get(gca,'yaxis'); set(gca,'YTickLabels',ya.TickValues*100)
            set(gca,'yminortick','on')
            ylabel('Fraction of active synapses')
            xlabel('AP #');
        end
    end
    subplot(s);
    if genotypes(g) == 148; title('wild-type 0.75Hz'); elseif genotypes(g) == 225; title('tomosynKO 0.75Hz');end
    ya = get(gca,'yaxis'); set(gca,'YTickLabels',ya.TickValues*100)
    set(gca,'yminortick','on')
    ylabel('Fraction of active synapses')
    xlabel('AP #');
    %camroll(270)
end


%%  Plot average 0.75Hz AP shape in active synapse (optional: per week)
stimStarts_075 = [600 615 630 645 660];
frontTail_075 = 4; backTail_075 = 9;
colors = ['k','b']; genotypes = [148, 225];
for w = 1%[22 25 26] % use the array with week labels to plot for each week separately
    f1=figure('name',['week ' num2str(w) ' amplitudes 0.75Hz']); set(f1,'windowstyle','docked'); clf;
    f2=figure('name',['week ' num2str(w) ' shapes 0.75Hz']); set(f2,'windowstyle','docked'); clf;
    f3=figure('name',['week ' num2str(w) ' COVs 0.75Hz']); set(f3,'windowstyle','docked'); clf;
    for g = 1:length(genotypes)
        shift_x = g-1; shift_jitter = 10000; color = colors(g);%clf;
        thisCondition = tomVSwt(tomVSwt.Condition == num2str(genotypes(g)),:);% & tomVSwt.Week == w,:); % use the extra filter option to plot for each week separately
        for ap = 1:length(stimStarts_075)
            % clear variables here if you plot them for per AP for all cells
            ampMeans = []; covMeans = [];  decayTimes = []; riseTimes = []; allShapesMeans = [];
            for cellNum = 1:size(thisCondition,1) 
                % clear variables here if you plot them per cell
                ampsMat = []; shapesMat = [];
                cellRaster = thisCondition.binaryRaster075{cellNum,1};
                activeSyn = find(cellRaster(:,ap) > 0);
                if ~isempty(activeSyn)
                    for idx = 1:length(activeSyn)
                        sig = thisCondition.dFoverF0_detrend{cellNum,1}(activeSyn(idx),stimStarts_075(ap)-frontTail_075:stimStarts_075(ap)+backTail_075);
                        shapesMat(end+1,:) = sig;
                        ampsMat(end+1,1) = thisCondition.fiveAt075Hz_dFoF0{cellNum,1}(activeSyn(idx),ap);
                    end
                
                    xpos = randi([(ap*10000)-2000 (ap*10000)+2000]+(shift_x*shift_jitter))/10000;
                    figure(f1);
                    subplot(1,5,ap); 
                    % plot 1 single mean amplitude of all active synapses for this AP/cell/week
                    plot(xpos,mean(ampsMat),['o' color]); hold on 
                    
                    figure(f3);
                    subplot(1,5,ap); 
                    % plot 1 single CoV of all active synapses for this AP/cell/week
                    plot(xpos,std(ampsMat)/mean(ampsMat),['o' color]); hold on 

                    %  store this mean for plotting overall median of all cells for this AP/week
                    ampMeans(end+1,1) = mean(ampsMat); 
                    %  store this COV for plotting overall median of all cells for this AP/week
                    covMeans(end+1,1) = std(ampsMat)/mean(ampsMat); 
                    try
                        % decay of average AP trace from active synapses from this cell/week
                        decayTimes(end+1,1) = -1/getApDecayRate(mean(shapesMat,1), 'ap','exp1st','shift'); 
                        % rise of average AP trace from active synapses from this cell/week
                        riseTimes(end+1,1) = -1/getApDecayRate(flip(mean(shapesMat,1)), 'ap','exp1st','shift'); 
                    catch ME
                        ...
                    end
                    allShapesMeans(end+1:end+size(shapesMat,1),:) = shapesMat;
                end
            end
            figure(f1);
            subplot(1,5,ap);
            meanBarAmp = mean(ampMeans(~isnan(ampMeans) & ampMeans ~= -Inf & ampMeans ~= Inf));
            plot([ap-0.4 ap+0.4]+shift_x, [meanBarAmp,meanBarAmp],['-' color])
            ylim([-0.05 1.2]); set(gca,'yminortick','on');
            xlim([ap-0.5 ap+0.5+1])
            
            figure(f3);
            subplot(1,5,ap);
            medianBarCov = median(covMeans(~isnan(covMeans) & covMeans ~= -Inf & covMeans ~= Inf));
            plot([ap-0.4 ap+0.4]+shift_x, [medianBarCov,medianBarCov],['-' color])
            ylim([-0.05 1]); 
            xlim([ap-0.5 ap+0.5+1])            

            figure(f2);
            subplot(3,5,ap);
            stdShade(allShapesMeans(:,3:end-2),0.1,[0 0 shift_x]); hold on;
            set(gca,'yminortick','on')
            ylim([-0.15 1.05]); xlim([0 12])

            subplot(3,5,[ap+5 ap+10]);
            plot(randi([(ap*10000)-3000 (ap*10000)-1000]+(shift_x*shift_jitter), [1, length(decayTimes)])/10000, decayTimes', ['o' color]); hold on
            plot([ap-0.35 ap-0.05]+shift_x, [median(decayTimes),median(decayTimes)],['-' color])
            ylabel('time (s)'); ylim([-0.1 3])
            set(gca,'yminortick','on')
            plot(randi([(ap*10000)+1000 (ap*10000)+3000]+(shift_x*shift_jitter), [1, length(riseTimes)])/10000, riseTimes',['d' color]); hold on  
            plot([ap+0.05 ap+0.35]+shift_x, [median(riseTimes),median(riseTimes)],['-' color])
            ylabel('time (s)'); ylim([-0.1 3])
            set(gca,'yminortick','on')
        end
    end
%      orient(f1,'landscape'); figure(f1); print(['allActiveAmps_' num2str(w) '_1'] ,'-fillpage','-dpdf','-painters');
%      orient(f2,'landscape'); figure(f2); print(['allActiveShapes_' num2str(w) '_2'],'-fillpage','-dpdf','-painters');
%      orient(f3,'landscape'); figure(f3); print(['allActiveCOVs_' num2str(w) '_3'],'-fillpage','-dpdf','-painters');
end

%%  Plot average 2Hz AP shape in active synapse (per week)
stimStart_2 = 810;
frontTail_2 = 3; backTail_2 = 6;
colors = ['k','b']; genotypes = [148, 225];
for w = 1%[22 25 26] % use the array with week labels to plot for each week separately
    f1=figure('name',['week ' num2str(w) ' 2Hz']); set(f1,'windowstyle','docked'); clf;
    for g = 1:length(genotypes)
        shift_x = g-1; shift_jitter = 10000; color = colors(g);%clf;
        thisCondition = tomVSwt(tomVSwt.Condition == num2str(genotypes(g)),:);% & tomVSwt.Week == w,:); % use the extra filter option to plot for each week separately
        ap=1;
        % clear variables here if you plot them for per AP for all cells
        ampMeans = []; covMeans = [];  decayTimes = []; riseTimes = []; allShapesMeans = [];
        for cellNum = 1:size(thisCondition,1) 
            % clear variables here if you plot them per cell
            ampsMat = []; shapesMat = [];
            cellRaster = thisCondition.binaryRaster2{cellNum,1};
            activeSyn = find(cellRaster(:,ap) > 0);
            if ~isempty(activeSyn)
                for idx = 1:length(activeSyn)
                    sig = thisCondition.dFoverF0_detrend{cellNum,1}(activeSyn(idx),stimStart_2-frontTail_2:stimStart_2+backTail_2);
                    shapesMat(end+1,:) = sig;
                    ampsMat(end+1,1) = thisCondition.fiveAt2Hz_dFoF0{cellNum,1}(activeSyn(idx),ap);
                end

                xpos = randi([(ap*10000)-2000 (ap*10000)+2000]+(shift_x*shift_jitter))/10000;
                s1=subplot(3,3,[1 4 7]); 
                % plot 1 single mean amplitude of all active synapses for this AP/cell/week
                plot(xpos,mean(ampsMat),['o' color]); hold on  
                s2=subplot(3,3,[2 5 8]); 
                % plot 1 single CoV of all active synapses for this AP/cell/week
                plot(xpos,std(ampsMat)/mean(ampsMat),['o' color]); hold on 

                %  store this mean for plotting overall median of all cells for this AP/week
                ampMeans(end+1,1) = mean(ampsMat); 
                %  store this COV for plotting overall median of all cells for this AP/week
                covMeans(end+1,1) = std(ampsMat)/mean(ampsMat); 
                try
                    % decay of average AP trace from active synapses from this cell/week
                    decayTimes(end+1,1) = -1/getApDecayRate(mean(shapesMat,1), 'ap','exp1st','shift'); 
                    % rise of average AP trace from active synapses from this cell/week
                    riseTimes(end+1,1) = -1/getApDecayRate(flip(mean(shapesMat,1)), 'ap','exp1st','shift'); %rise
                catch ME
                    ...
                end
                allShapesMeans(end+1:end+size(shapesMat,1),:) = shapesMat;
            end
        end
        subplot(s1);
        meanBarAmp = mean(ampMeans(~isnan(ampMeans) & ampMeans ~= -Inf & ampMeans ~= Inf));
        plot([ap-0.4 ap+0.4]+shift_x, [meanBarAmp,meanBarAmp],['-' color])
        ylim([-0.05 0.8]); set(gca,'yminortick','on');
        xlim([ap-0.5 ap+0.5+1])

        subplot(s2);
        medianBarCov = median(covMeans(~isnan(covMeans) & covMeans ~= -Inf & covMeans ~= Inf));
        plot([ap-0.4 ap+0.4]+shift_x, [medianBarCov,medianBarCov],['-' color])
        %ylim([-0.05 1]); 
        xlim([ap-0.5 ap+0.5+1])            

        s3=subplot(333);
        stdShade(allShapesMeans(:,3:end-2),0.1,[0 0 shift_x]); hold on;
        set(gca,'yminortick','on')
        ylim([-0.15 1.05]); xlim([0 8])

        s3=subplot(3,3,[6 9]);
        plot(randi([(ap*10000)-3000 (ap*10000)-1000]+(shift_x*shift_jitter), [1, length(decayTimes)])/10000, decayTimes', ['o' color]); hold on
        plot([ap-0.35 ap-0.05]+shift_x, [median(decayTimes),median(decayTimes)],['-' color])
        ylabel('time (s)'); ylim([-0.1 5.5])
        set(gca,'yminortick','on')
        plot(randi([(ap*10000)+1000 (ap*10000)+3000]+(shift_x*shift_jitter), [1, length(riseTimes)])/10000, riseTimes',['d' color]); hold on  
        plot([ap+0.05 ap+0.35]+shift_x, [median(riseTimes),median(riseTimes)],['-' color])
        ylabel('time (s)'); ylim([-0.1 5.5])
        set(gca,'yminortick','on')

    end
%     orient(f1,'landscape'); figure(f1); print(['allActiveAmps_' num2str(w) '_2Hz'] ,'-fillpage','-dpdf','-painters');
end



%% 2Hz  AP ratios: 2nd/1st and 5th/1st
colors = ['k','b']; genotypes = [148, 225]; clf; 
shift=0;
wt=mean(cell2mat(table2cell(tomVSwt(tomVSwt.Condition == '148',{'fiveAt2Hz_dFoF0'}))),1);  % mean value from all wt measurements,
%     potentially to be used for normalization
for ap = [2,5]
    apRatios = cellfun(@(x) x(:,ap)./x(:,1), tomVSwt.fiveAt2Hz_dFoF0, 'uniformoutput', false);
    for g = 1:length(genotypes)
        aprMeans = cell2mat(cellfun(@(x) mean(x), apRatios(tomVSwt.Condition == num2str(genotypes(g))), 'uniformoutput', false));
        aprMeans = aprMeans(~isnan(aprMeans) & aprMeans~= -Inf & aprMeans~= Inf);

        hold on;
        if ap==5; mType = 's'; else; mType = 'o'; end
        scatter(randi([-200 200]+(shift*1000), [1,length(aprMeans)])/1000, aprMeans, 72, [mType colors(g)]);
        plot([shift-0.3 shift+0.3], [median(aprMeans) median(aprMeans)], ['-' colors(g)])
        shift=shift+1;
        ylim([-1 4])
    end
end

%% Examples of individual traces
% to do (for final presentation)

%% Examples of raster plots for 0.75Hz
% to do (for final presentation)

%% RRP depletion during 40Hz stimulation
genotypes = [148, 225];
f=figure('name','RRP depletion rate'); set(f,'windowstyle','docked'); clf;
for g = 1:length(genotypes)
    allCellsMeans = cell2mat(cellfun(@(x) mean(x,1), table2cell(tomVSwt(tomVSwt.Condition == num2str(genotypes(g)),'dFoverF0_detrend')), 'uniformoutput', false));
    stimInt = tomVSwt.hundredAt40Hz_frameNum{1,1};
    for meanTrace = 1:size(allCellsMeans,1)
        train = allCellsMeans(meanTrace,:);
        trainPeak = max(train([1020 1021 1022 1024 1025]));

        trainPeakIdx = find(train==trainPeak,1);
        trainEndIdx = find(train==max(train([1045:1048])),1);

        trainDecay = train([trainPeakIdx:trainEndIdx]);
        [d, gof] = getApDecayRate(trainDecay,'40Hz','exp2nd','shift');
        decs_stim(meanTrace,[1 2]+(2*g-2)) = -1./d;
        goods_stim(meanTrace,1+g-1) = gof;

        probeDecay = train([trainEndIdx:1054]);
        [d, gof] = getApDecayRate(probeDecay,'40Hz','exp1st','shift');
        decs_probe(meanTrace,1+g-1) = -1./d;
        goods_probe(meanTrace,1+g-1) = gof;     
    end
end

jit1 = randi([920, 1070], [5,length(decs_stim(decs_stim(:,1)~=0,1))])/1000;
jit2 = randi([1920, 2070], [5,length(decs_stim(:,3))])/1000;

subplot(151)
boxplot([decs_stim(decs_stim(:,1)~=0,1);  decs_stim(:,3)],[ones([1, length(jit1)]), ones([1, length(jit2)])+1],'notch','on');
hold on; 
scatter(jit1(1,:), decs_stim(decs_stim(:,1)~=0,1),'ok');  
scatter(jit2(1,:), decs_stim(:,3),'ob'); 
set(gca,'xticklabel',{'wt','tomKO'})
title({'1st component of', 'the RRP depletion rate'})
ylabel('time (s)'); %ylim([-3 11])

subplot(152)
boxplot([decs_stim(decs_stim(:,2)~=0,2);  decs_stim(:,4)],[ones([1, length(jit1)]), ones([1, length(jit2)])+1],'notch','on');
hold on; 
scatter(jit1(2,:), decs_stim(decs_stim(:,2)~=0,2),'ok');  
scatter(jit2(2,:), decs_stim(:,4),'ob'); 
set(gca,'xticklabel',{'wt','tomKO'})
title({'2nd component of', 'the RRP depletion rate'})
ylabel('time (s)'); ylim([4 25])

subplot(153)
boxplot([goods_stim(goods_stim(:,1)~=0,1);  goods_stim(:,2)],[ones([1, length(jit1)]), ones([1, length(jit2)])+1],'notch','on');
hold on; 
scatter(jit1(3,:), goods_stim(goods_stim(:,1)~=0,1),'dk');  
scatter(jit2(3,:), goods_stim(:,2),'db'); 
set(gca,'xticklabel',{'wt','tomKO'})
title({'goodness of fit', 'for depletion rates'})
ylabel('time (s)');
xlim([0 3])

subplot(154)
boxplot([decs_probe(decs_probe(:,1)~=0,1);  decs_probe(:,2)],[ones([1, length(jit1)]), ones([1, length(jit2)])+1],'notch','on');
hold on; 
scatter(jit1(4,:), decs_probe(decs_probe(:,1)~=0,1),'ok');  
scatter(jit2(4,:), decs_probe(:,2),'ob'); 
set(gca,'xticklabel',{'wt','tomKO'})
title('probe decay rate')
ylabel('time (s)'); ylim([1 7]) % some outliers get excluded

subplot(155);hold on
scatter(jit1(5,:), goods_probe(decs_probe(:,1)~=0,1),'dk');  
scatter(jit2(5,:), goods_probe(:,2),'db'); 
set(gca,'xticklabel',{'','tomKO','wt',''})
title({'goodness of fit', 'for probe dissociation rates'})
ylabel('time (s)'); 
xlim([0 3])


%% RRP replenishment after 40Hz stimulation
% load already existing data;
load tomVswt_rasterRecovery.mat
tomVSwt = [tomVSwt, tomVswt_rasterRecovery];

% first detect active recovery potentials
% stimStarts = [1055 1085];
% frontTail = 2;
% backTail = 6;

%%% this code runs the detection from scratch
% parfor cellNum = 1:size(tomVSwt,1)
%     apScores = [];
%     apPositions = [];
%     apBinary = [];
%     for synNum = 1:size(tomVSwt.dFoverF0_detrend{cellNum,1},1)
%         for apTime = 1:length(stimStarts)
%             noiseTrace = tomVSwt.dFoverF0_detrend{cellNum,1}(synNum,1065:1075);
%             apTrace = tomVSwt.dFoverF0_detrend{cellNum,1}(synNum,stimStarts(apTime)-frontTail:stimStarts(apTime)+backTail);
%             
%             lagFrames = [frontTail+1:frontTail+3];
%             peakLag = cell2mat(findStimulationLag(apTrace(lagFrames)));
%             peakIdx = lagFrames(peakLag+length(lagFrames)-1);
%             
%             shiftFactor = min(apTrace);     % apTrace is being shifted above zero just for detection purposes
%             if shiftFactor > 0
%                 apTrace = apTrace-abs(shiftFactor);
%             elseif shiftFactor < 0
%                 apTrace = apTrace+abs(shiftFactor);
%             end
%             
%             [apBin, apScore, apPos] = isItAP(apTrace, noiseTrace, 'mini',5,peakIdx); 
%             
%             apScores(synNum, apTime) = apScore;
%             if apPos ~= 0
%                apPos = apPos + (stimStarts(apTime)-frontTail-1);
%             end 
%             apPositions(synNum, apTime) = apPos;
%             apBinary(synNum, apTime) = apBin;
%         end
%         
%     end
%     allApEvals_bin{cellNum,1} = apBinary;
%     allApEvals_scores{cellNum,1} = apScores;
%     allApEvals_pos{cellNum,1} = apPositions;
%     %return
% end
% fprintf('Done.\n')
% tomVSwt = [tomVSwt table(allApEvals_bin, allApEvals_scores, allApEvals_pos,'variablenames',{'binaryRasterRecovery','apScoresRecovery','apPositionsRecovery'})];

% then plot ratios of active recovery APs
clearvars -except tom*
genotypes = [148, 225];
f=figure('name','RRP replenishment rate'); set(f,'windowstyle','docked'); clf;
stimStarts = [1055 1085];
colors = ['k','b'];
for g = 1:length(genotypes)
    ratMeans=[];
    for ap = 1:length(stimStarts)
        thisCondition = tomVSwt(tomVSwt.Condition == num2str(genotypes(g)),:);
        r=1;
        for cellNum = 1:size(thisCondition,1)
            cellRaster = thisCondition.binaryRasterRecovery{cellNum,1};
            activeSyn = find(cellRaster(:,ap)==1); % go only through active synapses
%             activeSyn = 1:length(cellRaster); % uncomment to go throught all synapses
            promRatios=[];
            if ~isempty(activeSyn)
                for aSyn = 1:length(activeSyn)
                    apProminence = diff([min(thisCondition.dFoverF0_detrend{cellNum,1}(activeSyn(aSyn),[stimStarts(ap)-1:stimStarts(ap)+1])), ...
                        max(thisCondition.dFoverF0_detrend{cellNum,1}(activeSyn(aSyn),[stimStarts(ap)-1:stimStarts(ap)+1]))]);
                    
                    fortyProminence = diff([min(thisCondition.dFoverF0_detrend{cellNum,1}(activeSyn(aSyn),[1019:1023])), ...
                        max(thisCondition.dFoverF0_detrend{cellNum,1}(activeSyn(aSyn),[1019:1023]))]);
                    
                    if apProminence>0 && fortyProminence>0
                        promRatios(end+1,1) = apProminence/fortyProminence;
                    end
                end
                ratMeans(r,ap) = mean(promRatios);
                r=r+1;
            end
        end
    end
%     ratiosAp21 = ratMeans(ratMeans(:,1)~=0,1);
%     ratiosAp51 = ratMeans(ratMeans(:,2)~=0,2);
    ratiosAp21 = ratMeans(:,1);
    ratiosAp51 = ratMeans(:,2);
    
    subplot(121+g-1); hold on;
    boxplot([ratiosAp21; ratiosAp51], [ones(size(ratiosAp21)); ones(size(ratiosAp51))+1],'notch','on');
%     scatter(randi([920, 1070], [1,length(ratiosAp1)])/1000,ratiosAp1,'k');
%     scatter(randi([920, 1070]+1000, [1,length(ratiosAp2)])/1000,ratiosAp2,'b');
    plot([1.1 1.9], [ratiosAp21 ratiosAp51],['-o' colors(g)])
    ylim([-0.05 3]); set(gca,'yminortick','on')
end

%% mini detection (about 10 min)
tomVSwt(strcmp(tomVSwt.FileName,{'Results_210622_TomKOvsWT_cs04_001'}),:) = []; % this cell is corrupted in baseline
clearvars -except tom*
miniFrontTail = 4;
miniBackTail = 12;
warning('off')
tic
parfor cellNum = 1:size(tomVSwt,1)
    allPos = []; allBin = [];
    scoreTreshold = 10; minPeakH = 0.07;

    for synNum = 1:size(tomVSwt.dFoverF0_detrend{cellNum,1},1)
        [~, pks] = findpeaks(tomVSwt.dFoverF0_detrend{cellNum,1}(synNum,10:590),'minpeakheight',minPeakH,'minpeakdistance',10);
        if ~isempty(pks) 
            pks=pks+9;
            for thisPeak = 1:length(pks)
                miniNoise = tomVSwt.dFoverF0_detrend{cellNum,1}(synNum,5:20);
                [isItEvent, ~, miniPosition] = isItAP(tomVSwt.dFoverF0_detrend{cellNum,1}(synNum,...
                    pks(thisPeak)-miniFrontTail:pks(thisPeak)+miniBackTail),miniNoise, 'mini',scoreTreshold);
                if miniPosition ~= 0
                   miniPosition = miniPosition + (pks(thisPeak)-miniFrontTail-1);
                end 
                allPos(synNum, pks(thisPeak)) = miniPosition;
                allBin(synNum, pks(thisPeak)) = isItEvent;
            end
        end
    end
    allMiniBinary{cellNum,1} = allBin;
    allMiniPosition{cellNum,1} = allPos;
end
warning('on')
toc
fprintf('Done.\n')
tomVSwt = [tomVSwt table(allMiniBinary, allMiniPosition,'variablenames',{'binaryMinis','apPositionsMini'})];

%%% scanning approach (slow, hours)
% scanFrameStart = 1;
% backTail = 7;
% clf
% 
% parfor cellNum = 1:size(tomVSwt,1)
%     apScores = [];
%     apPositions = [];
%     apBinary = [];
%     for synNum = 1:size(tomVSwt.dFoverF0_detrend{cellNum,1},1)
% %         plot(tomVSwt.dFoverF0_detrend{cellNum,1}(synNum,1:800),'k'); hold on;
%         for scanFrameIdx = scanFrameStart:599-scanFrameStart-backTail
%             noiseStart = 1;%randi(589);
%             noiseTrace = tomVSwt.dFoverF0_detrend{cellNum,1}(synNum,noiseStart:noiseStart+backTail);
%             [apBin, apScore, apPos] = isItAP(tomVSwt.dFoverF0_detrend{cellNum,1}(synNum,scanFrameIdx:scanFrameIdx+backTail), noiseTrace, 'mini'); 
%             
% %             plot(scanFrameIdx:scanFrameIdx+backTail,tomVSwt.dFoverF0_detrend{cellNum,1}(synNum,scanFrameIdx:scanFrameIdx+backTail));
% %             if apBin > 0 && apPos > 0
% %                 plot(apPos+scanFrameIdx-1,tomVSwt.dFoverF0_detrend{cellNum,1}(synNum,apPos+scanFrameIdx-1),'ok');
% %             end
%             
%             apScores(synNum, scanFrameIdx+2) = apScore;
%             if apPos ~= 0
%                apPos = apPos + scanFrameIdx-2;
%             end 
%             apPositions(synNum, scanFrameIdx+2) = apPos;
%             apBinary(synNum, scanFrameIdx+2) = apBin;
%         end
%         
%     end
%     allApEvals_bin{cellNum,1} = apBinary;
%     allApEvals_scores{cellNum,1} = apScores;
%     allApEvals_pos{cellNum,1} = apPositions;
%     %return
% end
% fprintf('Done.\n')

%% mini quantification and visualisation
clearvars -except tom*
genotypes = [148, 225];

for g=1:length(genotypes)
    thisCondition = tomVSwt(tomVSwt.Condition == num2str(genotypes(g)),:);
    for cellNum = 1:size(thisCondition,1)
        cellRaster = thisCondition.binaryMinis{cellNum,1};
        [activeSyn, miniTime] = find(cellRaster==1);
        if ~isempty(activeSyn)
            m=1; miniShapes=[]; miniAmps=[];
            for miniNum = 1:length(activeSyn)
                mini = thisCondition.dFoverF0_detrend{cellNum,1}(activeSyn(miniNum),miniTime(miniNum)-3:miniTime(miniNum)+7);
                miniShapes(m,:) = mini;
                miniAmps(m,1) = thisCondition.dFoverF0_detrend{cellNum,1}(activeSyn(miniNum),miniTime(miniNum));
                m=m+1;
            end
        end
        allMiniShapes{cellNum,g} = miniShapes;
        allMiniAmps{cellNum,g} = miniAmps;
    end
    
    % get percentage of synapses that show spontaneous release 
    miniFreq = cellfun(@(x, y) (length(x)./size(y,1))*100,allMiniAmps(1:size(thisCondition,1),g),thisCondition.RawSignal,'uniformoutput',false);
    miniFreqs{1,g} = miniFreq;
end
f=figure('name','mini histograms'); set(f,'windowstyle','docked'); clf;
subplot(2,3,[1:3])
histogram(cell2mat(allMiniAmps(:,2)),[0:0.01:0.5]); hold on;
histogram(cell2mat(allMiniAmps(:,1)),[0:0.01:0.5])
ylabel('mGT count')
xlabel('mGT amplitude')

clearvars -except tom* *Mini* genotypes
ampsPerCell=[]; shapesPerCell = [];
for g = 1:length(genotypes)
    thisCondition = tomVSwt(tomVSwt.Condition == num2str(genotypes(g)),:);
    ap = 1; 
    for cellNum = 1:size(thisCondition,1) 
        stims = {'075';'2'};
        for stim = 1:length(stims)
            if stim == 1; stim_start = 600; frontTail = 3; backTail = 10; elseif stim == 2; stim_start = 810; frontTail = 1; backTail = 4; end
            ampsMat = []; shapesMat = [];
            cellRaster = thisCondition.(['binaryRaster' stims{stim,1}]){cellNum,1};
            activeSyn = find(cellRaster(:,ap) > 0);
            if ~isempty(activeSyn)
                for idx = 1:length(activeSyn)
                    ampsMat(end+1,1) = thisCondition.(['fiveAt' stims{stim,1} 'Hz_dFoF0']){cellNum,1}(activeSyn(idx),ap);
                    sig = thisCondition.dFoverF0_detrend{cellNum,1}(activeSyn(idx),stim_start-frontTail:stim_start+backTail);
                    shapesMat(end+1,:) = sig;
                end
            end
            ampsPerCell(cellNum,stim+(2*(g-1))) = mean(ampsMat);
            shapesPerCell{cellNum,stim+(2*(g-1))} = shapesMat;
        end
    end
    subplot(2,3,3+stim-1); hold on % 0.75 Hz
    stdShade(cell2mat(shapesPerCell(:,end-1)),0.1,[0 0 g-1]); hold on;
    set(gca,'yminortick','on')
    ylim([-0.15 0.9]);
    legend({'wt std','wt mean','tomKO std','tomKO mean'})
    title('1st evoked AP in 0.75 Hz')
    xlabel('frames')
    ylabel('iGluSnFR amplitude')
    
    subplot(2,3,3+stim); hold on % 2 Hz
    stdShade(cell2mat(shapesPerCell(:,end)),0.1,[0 0 g-1]); hold on;
    set(gca,'yminortick','on')
    ylim([-0.15 0.9]); xlim([0 7])
    legend({'wt std','wt mean','tomKO std','tomKO mean'})
    title('1st evoked AP in 2 Hz')
    xlabel('frames')
    ylabel('iGluSnFR amplitude')
end

subplot(2,3,6); hold on % 2 Hz
stdShade(cell2mat(allMiniShapes(:,1)),0.1,[0 0 0]); hold on;
stdShade(cell2mat(allMiniShapes(:,2)),0.1,[0 0 1]);
set(gca,'yminortick','on')
ylim([-0.15 0.9]); 
legend({'wt std','wt mean','tomKO std','tomKO mean'})
title('average mini shape')
xlabel('frames')
ylabel('iGluSnFR amplitude')

mean075wt=ampsPerCell(~isnan(ampsPerCell(:,1)) & ampsPerCell(:,1) ~= 0,1);
mean075tom=ampsPerCell(~isnan(ampsPerCell(:,3)) & ampsPerCell(:,3) ~= 0,3);
mean2wt=ampsPerCell(~isnan(ampsPerCell(:,2)) & ampsPerCell(:,2) ~= 0,2);
mean2tom=ampsPerCell(~isnan(ampsPerCell(:,4)) & ampsPerCell(:,4) ~= 0,4);

subplot(2,3,[1:3])
plot([mean(mean075wt) mean(mean075wt)],[0 200], 'k');
plot([mean(mean2wt) mean(mean2wt)],[0 200], '-.k');
plot([mean(mean075tom) mean(mean075tom)],[0 200], 'b');
plot([mean(mean2tom) mean(mean2tom)],[0 200], '-.b');
legend({'tomKO','wt','wt 0.75Hz','wt 2Hz','tomKO 0.75Hz', 'tomKO 2Hz'})


