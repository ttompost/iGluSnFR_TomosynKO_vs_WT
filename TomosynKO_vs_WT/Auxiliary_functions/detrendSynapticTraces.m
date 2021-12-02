function out  = detrendSynapticTraces(tracesTable, expType, varargin)
%%% The function detrends raw deltaF/F0 synaptic traces according to the input window
% sizes. Length of recordings is appropriated based on the experiment 
% type (different ca2+ concentrations or tomosyn versus wild-type).

%%% Author: Tea Tompos (teatompos@gmail.com)
% Last edited: 21st august 2021

    switch expType % number of frames per recording
        case 'tomVSwt'
            rec1 = 1:800;   
            rec2 = 801:1010;   
            rec3 = 1011:1110; 
        case 'diffBathCa'
            rec1 = 1:1100;      
            rec2 = 1101:1200;  
            rec3 = 1201:1300;  
        otherwise
            disp('unknown experiment type')
            return
    end
    varargin = varargin{1, 1};
    
    switch size(varargin,2)
        case 0 % do not detrend anything
            trace1 = cellfun(@(x) x(:,rec1), tracesTable.dFoverF0,'uniformoutput',false); 
            trace2 = cellfun(@(x) x(:,rec2), tracesTable.dFoverF0,'uniformoutput',false); 
            trace3 = cellfun(@(x) x(:,rec3), tracesTable.dFoverF0,'uniformoutput',false); 
        case 1
            winRec1 = varargin(1);
            trace1 = cellfun(@(x) x(:,rec1) - movmedian(x(:,rec1),winRec1,2), tracesTable.dFoverF0,'uniformoutput',false); 
            trace2 = cellfun(@(x) x(:,rec2), tracesTable.dFoverF0,'uniformoutput',false); 
            trace3 = cellfun(@(x) x(:,rec3), tracesTable.dFoverF0,'uniformoutput',false); 
        case 2
            winRec1 = varargin(1); winRec2 = varargin(2);
            trace1 = cellfun(@(x) x(:,rec1) - movmedian(x(:,rec1),winRec1,2), tracesTable.dFoverF0,'uniformoutput',false); 
            trace2 = cellfun(@(x) x(:,rec2) - movmedian(x(:,rec2),winRec2,2), tracesTable.dFoverF0,'uniformoutput',false);
            trace3 = cellfun(@(x) x(:,rec3), tracesTable.dFoverF0,'uniformoutput',false); 
        case 3
            winRec1 = varargin(1); winRec2 = varargin(2); winRec3 = varargin(3);
            trace1 = cellfun(@(x) x(:,rec1) - movmedian(x(:,rec1),winRec1,2), tracesTable.dFoverF0,'uniformoutput',false); 
            trace2 = cellfun(@(x) x(:,rec2) - movmedian(x(:,rec2),winRec2,2), tracesTable.dFoverF0,'uniformoutput',false); 
            trace3 = cellfun(@(x) x(:,rec3) - movmedian(x(:,rec3),winRec3,2), tracesTable.dFoverF0,'uniformoutput',false);
    end

    traces = cellfun(@(x,y,z) horzcat(x,y,z), trace1, trace2, trace3, 'uniformoutput',false);
    out = [tracesTable, cell2table(traces,'VariableNames',{'dFoverF0_detrend'})];
end