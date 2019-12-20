
function configuration = SimulatedAnneling(im, configuration, f, miuf, alpha, beta, mask, T)
    % potential
    localPotential = calclocalPotentials(im, configuration, f, miuf, alpha, beta);
    % potential over current temperature
    prob = exp(-localPotential./T);
    % normalization
    probNorm = sum(prob,3);
    prob(:,:,1) = prob(:,:,1) ./ probNorm;
    for i = 2:size(localPotential,3) 
        prob(:,:,i) = prob(:,:,i) ./ probNorm;
        prob(:,:,i) = prob(:,:,i) + prob(:,:,i-1);
    end
    % random probability
    randProb = rand(size(probNorm));
    id = cumsum(prob>randProb,3);
    id = id(:,:,end);
    id = size(localPotential,3) - id;
    id = id + 1;
    % find the matched label and assign
    configuration(mask) = f(id(mask));
    
    % do again for another side of mask
    localPotential = calclocalPotentials(im, configuration, f, miuf, alpha, beta);
    prob = exp(-localPotential./T);
    probNorm = sum(prob,3);
    prob(:,:,1) = prob(:,:,1) ./ probNorm;
    for i = 2:size(localPotential,3) 
        prob(:,:,i) = prob(:,:,i) ./ probNorm;
        prob(:,:,i) = prob(:,:,i) + prob(:,:,i-1);
    end
    randProb = rand(size(probNorm));
    id = cumsum(prob>randProb,3);
    id = id(:,:,end);
    id = size(localPotential,3) - id;
    id = id + 1;
    
    configuration(~mask) = f(id(~mask));
end