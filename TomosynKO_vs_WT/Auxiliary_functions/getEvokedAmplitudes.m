function out = getEvokedAmplitudes(tracesTable, expType)
%%% The function extracts precise fluorescence peak values from individual synaptic
%%% traces during electrical stimulation. Stimulation timepoints are appropriated based 
% on the experiment type (different ca2+ concentrations or tomosyn versus wild-type).

%%% Author: Tea Tompos (teatompos@gmail.com)
% Last edited: 23rd august 2021
    
    try % use detrended deltaF/F0 traces if the input table contains those
       dataToAnalyse = table2cell(tracesTable(:,'dFoverF0_detrend'));
    catch ME % use raw deltaF/F0 traces if the input table contains those
       dataToAnalyse = table2cell(tracesTable(:,'dFoverF0'));   
    end
    
    out = tracesTable;
    switch expType % number of frames per recording
        case 'tomVSwt'
            stim2=[810:5:830]; % 5 ap at 2 hz
            backTail2 = 1; 
            lags2Hz = cellfun(@(x) cell2mat([findStimulationLag(x(:,stim2(1):stim2(1)+backTail2)),...
                findStimulationLag(x(:,stim2(2):stim2(2)+backTail2)),findStimulationLag(x(:,stim2(3):stim2(3)+backTail2)),...
                findStimulationLag(x(:,stim2(4):stim2(4)+backTail2)),findStimulationLag(x(:,stim2(5):stim2(5)+backTail2))]),...
                dataToAnalyse, 'uniformoutput', false);
            for cl = 1:size(dataToAnalyse,1)
                fiveAt2Hz_dfof{cl,1} = cell2mat(arrayfun(@(syn) dataToAnalyse{cl,1}(syn,stim2+lags2Hz{cl,1}(syn,:)), 1:size(dataToAnalyse{cl,1},1), 'uniformoutput', false)');        
                fiveAt2Hz_frameNum{cl,1} = lags2Hz{cl,1} + stim2;
            end
            out = [out, table(fiveAt2Hz_frameNum,fiveAt2Hz_dfof,'variablenames', {'fiveAt2Hz_frameNum','fiveAt2Hz_dFoF0'})];
            
            stim40 = [1021:1:1045]; % 100 ap at 40 hz 
            r1 = [1054:1:1056]; % 1 ap (0.2s)
            r2 = [1084:1:1086]; % 1 ap (0.2s)
        case 'diffBathCa'
            stim5 = [980:2:988];  % 5 ap at 5 hz (200ms each AP, or 2 frames for each AP)
            backTail5 = 1; 
            lags5Hz = cellfun(@(x) cell2mat([findStimulationLag(x(:,stim5(1):stim5(1)+backTail5)),...
                findStimulationLag(x(:,stim5(2):stim5(2)+backTail5)),findStimulationLag(x(:,stim5(3):stim5(3)+backTail5)),...
                findStimulationLag(x(:,stim5(4):stim5(4)+backTail5)),findStimulationLag(x(:,stim5(5):stim5(5)+backTail5))]),...
                dataToAnalyse, 'uniformoutput', false);
            for cl = 1:size(dataToAnalyse,1)
                fiveAt5Hz_dfof{cl,1} = cell2mat(arrayfun(@(syn) dataToAnalyse{cl,1}(syn,stim5+lags5Hz{cl,1}(syn,:)), 1:size(dataToAnalyse{cl,1},1), 'uniformoutput', false)');        
                fiveAt5Hz_frameNum{cl,1} = lags5Hz{cl,1} + stim5;
            end
            out = [out, table(fiveAt5Hz_frameNum,fiveAt5Hz_dfof,'variablenames', {'fiveAt5Hz_frameNum','fiveAt5Hz_dFoF0'})];
            
            stim20 = [1111:1:1114]; % 5 ap at 20 hz (300ms, or 3 frames for full simulation), actually 5 frames
            fiveAt20Hz_dfof = cellfun(@(x) x(:,stim20), dataToAnalyse, 'uniformoutput', false);
            fiveAt20Hz_frameNum = repelem(mat2cell(stim20,1,length(stim20)),size(dataToAnalyse,1),1); 
            out = [out, table(fiveAt20Hz_frameNum,fiveAt20Hz_dfof,'variablenames', {'fiveAt20Hz_frameNum','fiveAt20Hz_dFoF0'})];
            
            stim40 = [1211:1:1220]; % 40 ap at 40 hz (1s, or 10 frames for full), actually 12 frames
            r1 = [1229:1:1235]; % 1 ap (0.7s)
            r2 = [1259:1:1265]; % 1 ap (0.7s)
        otherwise
            disp('unknown experiment type')
            return
    end
    
    stim075=[600:15:660]; % 5 ap at 0.75 hz
    frontTail075 = 1; backTail075 = 1;
    % center peaks for 0.75Hz stimulation
    lags075Hz = cellfun(@(x) cell2mat([findStimulationLag(x(:,stim075(1)-frontTail075:stim075(1)+backTail075)),...
        findStimulationLag(x(:,stim075(2)-frontTail075:stim075(2)+backTail075)),...
        findStimulationLag(x(:,stim075(3)-frontTail075:stim075(3)+backTail075)),...
        findStimulationLag(x(:,stim075(4)-frontTail075:stim075(4)+backTail075)),...
        findStimulationLag(x(:,stim075(5)-frontTail075:stim075(5)+backTail075))]),...
        dataToAnalyse, 'uniformoutput', false);

    % center peaks for 1st recovery pulse 
    lags1ap = cellfun(@(x) cell2mat(findStimulationLag(x(:,r1))),dataToAnalyse, 'uniformoutput', false);
    
     % center peaks for 2nd recovery pulse 
    lags2ap = cellfun(@(x) cell2mat(findStimulationLag(x(:,r2))),dataToAnalyse, 'uniformoutput', false);
    
    for cl = 1:size(dataToAnalyse,1)
        fiveAt075Hz_dfof{cl,1} = cell2mat(arrayfun(@(syn) dataToAnalyse{cl,1}(syn,stim075+lags075Hz{cl,1}(syn,:)), 1:size(dataToAnalyse{cl,1},1), 'uniformoutput', false)');
        fiveAt075Hz_frameNum{cl,1} = lags075Hz{cl,1} + stim075;
        
        ap1_dfof{cl,1} = cell2mat(arrayfun(@(syn) dataToAnalyse{cl,1}(syn,r1(lags1ap{cl,1}(syn)+length(r1)-1)), 1:size(dataToAnalyse{cl,1},1), 'uniformoutput', false)'); 
        ap1_frameNum{cl,1} = r1(lags1ap{cl,1}+length(r1)-1)';
        
        ap2_dfof{cl,1} = cell2mat(arrayfun(@(syn) dataToAnalyse{cl,1}(syn,r2(lags2ap{cl,1}(syn)+length(r2)-1)), 1:size(dataToAnalyse{cl,1},1), 'uniformoutput', false)');
        ap2_frameNum{cl,1} = r2(lags1ap{cl,1}+length(r2)-1)';
    end
    
    hundredAt40Hz_dfof = cellfun(@(x) x(:,stim40), dataToAnalyse, 'uniformoutput', false);
    hundredAt40Hz_frameNum = repelem(mat2cell(stim40,1,length(stim40)),size(dataToAnalyse,1),1);

    out = [out, table(fiveAt075Hz_frameNum, fiveAt075Hz_dfof, hundredAt40Hz_frameNum, hundredAt40Hz_dfof,  ap1_frameNum, ap1_dfof,  ap2_frameNum, ap2_dfof,...
        'VariableNames', {'fiveAt075Hz_frameNum','fiveAt075Hz_dFoF0','hundredAt40Hz_frameNum','hundredAt40Hz_dFoF0',...
        'ap1_frameNum','ap1_dFoF0', 'ap2_frameNum','ap2_dFoF0'})];
end