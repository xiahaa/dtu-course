function configuration = ICM(im, configuration, f, miuf, alpha, beta, mask)
	% compute potentials    
    localPotential = calclocalPotentials(im, configuration, f, miuf, alpha, beta);
    % find the label which produces the minimum cost
    [~,id] = min(localPotential,[],3);
%     newconfiguration = configuration;
    configuration(mask) = f(id(mask));
    % do again for a new mask
    localPotential = calclocalPotentials(im, configuration, f, miuf, alpha, beta);
    [~,id] = min(localPotential,[],3);
    configuration(~mask) = f(id(~mask));
end
