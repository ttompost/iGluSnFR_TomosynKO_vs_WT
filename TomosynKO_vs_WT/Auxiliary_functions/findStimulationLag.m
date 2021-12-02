function stimLag = findStimulationLag(trace)
    for rw = 1:size(trace,1)
        maxIdx = find(max(trace(rw,:))==trace(rw,:),1);
        stimLag{rw,1} = maxIdx-(length(trace(rw,:))-1);
    end
end