function [eventScore, score, peakPos] = isItAP(trace, noise, type, varargin)
    % five step evaluation of a suspected action potential or miniature
    % glutamate transient (mGT) event
    
    % output: [apBinaryScore, apTotalScore, apPosition]
    
    warning('off')
    if length(trace) <= 5 
        disp('Trace is too short.. \n');
        return
    end
    v=cell2mat(varargin);
    if length(varargin)==1
        scoreThreshold = v(1);
    else
        scoreThreshold = 7;
    end
    
    if length(varargin)==2
        scoreThreshold = v(1);
        peakIdx = v(2);
        peakAmp = trace(peakIdx);
    else
        switch type 
            case 'mini'
                peakIdx = 3;
                peakAmp = trace(peakIdx);
            case 'ap'
                peakAmp = max(trace);
                peakIdx = find(trace==peakAmp,1);
        end
    end
    
    if length(noise) <= 4 
        disp('Noise trace is too short. \n');
        return
    end
    
    
    peakSnr = snr(peakAmp, mean(noise));
    traceGrad = gradient(trace);
    try peakGrad = traceGrad(peakIdx-1);
    catch
        peakGrad = [];
    end
    
    % evaluate detected APs
    score = 0;
    
    % score amplitude
    if  peakAmp >= 0.8
        score = score + 5; 
    elseif peakAmp < 0.8 && peakAmp >= 0.6
        score = score + 4;
    elseif peakAmp < 0.6 && peakAmp >= 0.4
        score = score + 3;
    elseif peakAmp < 0.4 && peakAmp >= 0.1
        score = score + 2;
    elseif peakAmp < 0.1 && peakAmp >= 0.07
        score = score + 1;
    elseif peakAmp < 0.07 
        score = score - 2;
    end
        
    % score snr
    if  peakSnr >= 50
        score = score + 3; 
    elseif peakSnr < 50 && peakSnr >= 30
        score = score + 2;
    elseif peakSnr < 30 && peakSnr >= 15
        score = score + 1;
    elseif peakSnr < 15 
        score = score - 1;
    end

    % score gradient
    if ~isempty(peakGrad)
        if peakGrad > 0.07
            score = score + 2; 
        elseif peakGrad < 0.07 && peakGrad >= 0.05
            score = score + 1; 
        end
    end
    
    % get shape points
    try peak_plus1 = trace(peakIdx+1);
    catch 
       peak_plus1 = [];
       score = score - 2;
    end
    
    try peak_plus2 = trace(peakIdx+2);
    catch 
       peak_plus2 = [];
       score = score - 2;
    end
    
    try peak_minus1 = trace(peakIdx-1);
    catch 
        peak_minus1 = [];
        score = score - 2;
    end
    
    try peak_minus2 = trace(peakIdx-2);
    catch
        peak_minus2 = [];
    end
    
    % get decay values
    decayPart = trace(peakIdx:end);
    if length(decayPart) >= 5
        [fitData, gofData] = fit([1:length(trace)]', trace', 'exp1', 'exclude', [1:peakIdx-1], 'startpoint', [1/mean(decayPart) -0.7]);
        coeffs = coeffvalues(fitData);
        decayRate = round(coeffs(2),3) * -1;
        goodFit = gofData.rsquare; %  goodnes of fit
        
        % score decay
        if 0.5 <= decayRate
            score = score + 1; 
        end 
        if goodFit >= 0.75
            score = score + 2;
        end
    else
        decayRate = [];
        goodFit = [];
        score = score - 2;
    end

    switch type
        case 'mini'
            if  ~isempty(peak_plus1) && ~isempty(peak_plus2) && ~isempty(peak_minus1)
                if peak_plus1 > peak_minus1
                    score = score + 1;
                else
                    score = score - 2;
                end
                if peak_plus1 > peak_plus2
                    score = score + 1;
                end      
                if peak_plus2 > peak_minus1
                    score = score + 1;
                end
                if peak_minus1 > peakAmp
                    score = score - 100;
                end                
            end
            
            % judge the event by its score
            if score >= scoreThreshold % max is 19, min is -11
                eventScore = 1;
            else
                eventScore = 0;
            end
            
        case 'ap'
            % score shape
            if  ~isempty(peak_plus1) && ~isempty(peak_plus2) && ~isempty(peak_minus1)
                if (peak_plus1 / peakAmp) < 0.35
                    if peak_minus1 > peak_plus1 
                        score = score + 1;
                    end
                    if ~isempty(peak_minus2) && peak_plus1 > peak_minus2
                        score = score + 1;
                    end 
                    if peak_minus1 > peak_plus2 
                        score = score + 1;
                    end                   
                elseif (peak_plus1 / peakAmp) >= 0.35
                    if peak_minus1 < peak_plus1 
                        score = score + 1;
                    end
                    if peak_plus1 > peak_plus2
                        score = score + 1;
                    end 
                    if ~isempty(peak_minus2) && peak_plus1 > peak_minus2 
                        score = score + 1;
                    end                     
                end
            end
            
            % judge the event by its score
            if score >= scoreThreshold % max is 19, min is -11
                eventScore = 1;
            else
                eventScore = 0;
            end
    end
    warning('on')

    if eventScore == 1
        peakPos = peakIdx;
    else
        peakPos = 0;
    end
    
end