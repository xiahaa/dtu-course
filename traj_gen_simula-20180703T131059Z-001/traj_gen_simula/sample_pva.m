function [pts,vts,ats,tss]=sample_pva(polyCoeffs, t, order)
    pts = [];
    vts = [];
    ats = [];
    numCoeff = order+1;
    tss = [];
    for i = 1:size(t,2)-1
        scalar = 1./(t(i+1)-t(i));
        ts = linspace(0,1,1000)';
        polyx = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,1);
        polyy = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,2);
        xsamples = polyx(1) + polyx(2).*ts + polyx(3).*ts.^2 + polyx(4).*ts.^3 + polyx(5).*ts.^4 + polyx(6).*ts.^5;
        ysamples = polyy(1) + polyy(2).*ts + polyy(3).*ts.^2 + polyy(4).*ts.^3 + polyy(5).*ts.^4 + polyy(6).*ts.^5;
        
        vxsamples = polyx(2) + polyx(3).*ts.*2 + polyx(4).*3.*ts.^2 + polyx(5).*4.*ts.^3 + polyx(6).*5.*ts.^4;
        vysamples = polyy(2) + polyy(3).*ts.*2 + polyy(4).*3.*ts.^2 + polyy(5).*4.*ts.^3 + polyy(6).*5.*ts.^4;
        
        axsamples = polyx(3).*2 + polyx(4).*6.*ts.^1 + polyx(5).*12.*ts.^2 + polyx(6).*20.*ts.^3;
        aysamples = polyy(3).*2 + polyy(4).*6.*ts.^1 + polyy(5).*12.*ts.^2 + polyy(6).*20.*ts.^3;
        
        vxsamples = vxsamples.*scalar;
        vysamples = vysamples.*scalar;
        axsamples = axsamples.*(scalar^2);
        aysamples = aysamples.*(scalar^2);
        
        pts = [pts;[xsamples ysamples]];
        vts = [vts;[vxsamples vysamples]];
        ats = [ats;[axsamples aysamples]];
        
        tss = [tss;ts+(i-1)];
    end
end