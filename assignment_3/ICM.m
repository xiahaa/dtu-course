function configuration = ICM(im, configuration, f, miuf, alpha, beta)
    mask = checkerboard(1,round(size(im,1)*0.5),round(size(im,2)*0.5));
    mask(size(im,1)+1:end,:) = [];
    mask(:,size(im,2)+1:end) = [];
    mask(mask > 0.5) = 1; mask = mask == 1;
    
    localPotential = calclocalPotentials(im, configuration, f, miuf, alpha, beta);
    [~,id] = min(localPotential,[],3);
%     newconfiguration = configuration;
    configuration(mask) = f(id(mask));
    
    localPotential = calclocalPotentials(im, configuration, f, miuf, alpha, beta);
    [~,id] = min(localPotential,[],3);
    configuration(~mask) = f(id(~mask));
end
