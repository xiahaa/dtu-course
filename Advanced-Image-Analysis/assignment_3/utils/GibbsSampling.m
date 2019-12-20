
function configuration = GibbsSampling(im, configuration, f, miuf, alpha, beta, mask)    
    newconfiguration = configuration;
    % firstly, compute potentials
    localPotential = calclocalPotentials(im, newconfiguration, f, miuf, alpha, beta);
    % potentials to probability
    prob = exp(-localPotential);
    % normalization
    probNorm = sum(prob,3);
    prob(:,:,1) = prob(:,:,1) ./ probNorm;
    for i = 2:size(localPotential,3) 
        prob(:,:,i) = prob(:,:,i) ./ probNorm;
        prob(:,:,i) = prob(:,:,i) + prob(:,:,i-1);
    end
    % random a probability
    randProb = rand(size(probNorm));
    % find the region where the random probability located inside
    id = cumsum(prob>randProb,3);
    id = id(:,:,end);
    id = size(localPotential,3) - id;
    id = id + 1;
    % assign the newconfiguration
    newconfiguration(mask) = f(id(mask));
    % do the previous again for another mask
    localPotential = calclocalPotentials(im, newconfiguration, f, miuf, alpha, beta);
    prob = exp(-localPotential);
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
    
    newconfiguration(~mask) = f(id(~mask));
    % if new cost is better, then take
    newCost = calcLikelihood(miuf,im,newconfiguration,alpha) + calcSmoothnessPrior(newconfiguration, beta);
    oldCost = calcLikelihood(miuf,im,configuration,alpha) + calcSmoothnessPrior(configuration, beta);
    
    if newCost < oldCost
        configuration = newconfiguration;
    end
end