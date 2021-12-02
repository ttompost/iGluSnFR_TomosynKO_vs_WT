function [decayRate, goodFit] = getApDecayRate(apTrace, stimType, fitType, normType)
    
    switch stimType
        case '40Hz'
            apPeakIdx = 1; % this is set as 1 because i pass only the decay part 
        case {'ap','mgt'}
            apPeakIdx = find(apTrace==max(apTrace),1); % used for 0.75Hz but its actually generic
        otherwise
            return
    end
    
    decayPart = apTrace(apPeakIdx:end);
    xData = (1:length(decayPart))+apPeakIdx-1;
    
    switch normType
        case 'shift'
            shiftFactor = min(decayPart);     % apTrace is being shifted above zero because fitnlm gives an error if there are negative values
            if shiftFactor > 0
                decayPart = decayPart-abs(shiftFactor);
            elseif shiftFactor < 0
                decayPart = decayPart+abs(shiftFactor);
            end
        case 'norm'
            normFactor = max(decayPart);  % normalization
            decayPart = decayPart/normFactor;
        case 'raw'
            % use raw trace  
            ...
    end
    
    switch fitType
        case 'exp1st'
            [fitData, gofData] = fit(xData', decayPart', 'exp1','startpoint', [1/mean(decayPart) -0.95]);
            coeffs = coeffvalues(fitData);
            decayRate = round(coeffs(2),3);
            
            %plot(fitData,xData, decayPart, 'predfunc'); 
            
            goodFit = gofData.rsquare;
            
        case 'exp2nd'
%             disp('unavailable at the moment')
%             return
            [fitData, gofData] = fit(xData', decayPart', 'exp2', 'startpoint', [1/mean(decayPart(1:2)) -0.7 1/mean(decayPart(3:end)) -0.25]); % 148: (3.7757   -0.556    3.618   -0.279) %% 225: (3662128 -0.886 1.224 -0.211)
            coeffs = coeffvalues(fitData);
            decayRate = [round(coeffs(2),3) round(coeffs(4),3)];
            
            %plot(fitData,xData, decayPart, 'predfunc'); 
            
            goodFit = gofData.rsquare;            
            
        case 'nlm' % generic non-linear model
            disp('unavailable at the moment')
            return
%             dataTable = table(xData(:), decayPart(:));
%             modelFunc = @(b,x) b(1) * exp(-b(2)*x(:, 1));  % Define the model as Y = a * exp(-b*x)
% 
%             aGuessed = 4; bGuessed = 0.1; % Guess values to start with.  Just make your best guess.
%             beta0 = [aGuessed, bGuessed];
% 
%             try 
%                 modParams = fitnlm(dataTable, modelFunc, beta0); % Compute the model and determine model parameters.
%                 coeffs = modParams.Coefficients{:, 'Estimate'}; % Extract the coefficient values from the the model object.
%                 yFitted = coeffs(1) * exp(-coeffs(2)*xData); % Create smoothed/regressed data using the model:
%                 decayRate = round(coeffs(2),3);
%             catch ME
%                 decayRate = 0;
%                 yFitted = zeros(size(decayPart));
%             end
% 
%             plot(apTrace, 'k'); hold on;
%             plot(xData+(apPeakIdx-1),decayPart,'ok');
%             plot(xData+(apPeakIdx-1),yFitted,'b','linewidth',1.6);
%             hold off
    end
end