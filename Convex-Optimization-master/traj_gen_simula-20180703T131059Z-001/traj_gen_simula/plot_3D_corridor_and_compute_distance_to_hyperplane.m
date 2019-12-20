function [h1, ts, ds1, ds2] = plot_3D_corridor_and_compute_distance_to_hyperplane(p, pz, l, polyCoeffs, realt)
    dims = 3;
    N = size(p,1);
    pWpts = zeros(N,dims);

    pWpts(:,1:2) = p;
    pWpts(:,3) = pz;
    
    %% compute hyperplane
    hyperplanes = zeros(size(pWpts,1)-1,8);
    
    supporting_points = [];
    
    for i = 1:1: size(pWpts,1)-1
        p31 = pWpts(i,:);
        p32 = pWpts(i+1,:);
        %% dir
        pdir = p32-p31;
        pdir = pdir./norm(pdir);
        %% n1
        n1 = [0 -pdir(3) pdir(2)];
        n1 = n1./norm(n1);
        d1 = -(n1(1)*pWpts(i,1)+n1(2)*pWpts(i,2)+n1(3)*pWpts(i,3));
        %% n2
        n2 = cross(n1,pdir);
        n2 = n2./norm(n2);
        d2 = -(n2(1)*pWpts(i,1)+n2(2)*pWpts(i,2)+n2(3)*pWpts(i,3));
        
        hyperplanes(i,:) = [n1 d1 n2 d2];
        
        length = norm(p32-p31);
        num = 10;
        samples = linspace(0,length,num);
        samples = samples';
        
        ptsSamples1 = repmat(p31,num,1) + samples.*repmat(pdir,num,1) - l.*repmat(n1,num,1) - l.*repmat(n2,num,1);
        ptsSamples2 = repmat(p31,num,1) + samples.*repmat(pdir,num,1) - l.*repmat(n1,num,1) + l.*repmat(n2,num,1);
        ptsSamples3 = repmat(p31,num,1) + samples.*repmat(pdir,num,1) + l.*repmat(n1,num,1) + l.*repmat(n2,num,1);
        ptsSamples4 = repmat(p31,num,1) + samples.*repmat(pdir,num,1) + l.*repmat(n1,num,1) - l.*repmat(n2,num,1);

        for j = 1:num
            h1 = plot3([ptsSamples1(j,1) ptsSamples2(j,1)],[ptsSamples1(j,2) ptsSamples2(j,2)],[ptsSamples1(j,3) ptsSamples2(j,3)],'--g','LineWidth', 2);
            h2 = plot3([ptsSamples2(j,1) ptsSamples3(j,1)],[ptsSamples2(j,2) ptsSamples3(j,2)],[ptsSamples2(j,3) ptsSamples3(j,3)],'--g','LineWidth', 2);
            h3 = plot3([ptsSamples3(j,1) ptsSamples4(j,1)],[ptsSamples3(j,2) ptsSamples4(j,2)],[ptsSamples3(j,3) ptsSamples4(j,3)],'--g','LineWidth', 2);
            h4 = plot3([ptsSamples4(j,1) ptsSamples1(j,1)],[ptsSamples4(j,2) ptsSamples1(j,2)],[ptsSamples4(j,3) ptsSamples1(j,3)],'--g','LineWidth', 2);
            h1.Color(4) = 0.5;
            h2.Color(4) = 0.5;
            h3.Color(4) = 0.5;
            h4.Color(4) = 0.5;
        end
    end
    
    segments = size(pWpts,1)-1;
    order = 5;
    numCoeff = order + 1;

    ts = [];
    ds1 = [];
    ds2 = [];
    
    for i=1:segments
        scalar = 1./(realt(i+1)-realt(i));
            
        polyx = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,1);
        polyy = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,2);
        polyz = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,3);
            
        hyperplane = hyperplanes(i,:);
        dist1 = hyperplane(1).*polyx + hyperplane(2).*polyy + hyperplane(3).*polyz;
        dist1(1) = dist1(1) + hyperplane(4);
        dist2 = hyperplane(5).*polyx + hyperplane(6).*polyy + hyperplane(7).*polyz;
        dist2(1) = dist2(1) + hyperplane(8);
        dist1 = flipud(dist1);
        dist2 = flipud(dist2);
        
        t1 = linspace(0,1,100)';
        ts = [ts;linspace(realt(i),realt(i+1),100)'];
                          
        for j = 1:numel(t1)
            d1 = polyval(dist1,t1(j));
            d2 = polyval(dist2,t1(j));
            ds1 = [ds1;abs(d1)];
            ds2 = [ds2;abs(d2)];
        end       
    end
end