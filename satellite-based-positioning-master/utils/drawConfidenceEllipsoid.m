function h = drawConfidenceEllipsoid(Dx, type, metrics)
%% modified code.
%% source code from aa@imm.dtu.dk, www.imm.dtu.dk/?aa
%% modified by xiahaa@space.dtu.dk
    % semi-axes in confidence ellipsoid
    % 95% fractile for 3 dfs is 7.815 = 2.796^2
    % 99% fractile for 3 dfs is 11.342 = 3.368^2
    %[vDxhat dDxhat] = eigsort(Dxhat(1:2,1:2));
    [vDxhat, dDxhat] = eigsort(Dx);
    %semiaxes = sqrt(diag(dDxhat));
    % 95% fractile for 2 dfs is 5.991 = 2.448^2
    % 99% fractile for 2 dfs is 9.210 = 3.035^2
    % df F(3,df).95 F(3,df).99
    % 1 215.71 5403.1
    % 2 19.164 99.166
    % 3 9.277 29.456
    % 4 6.591 16.694
    % 5 5.409 12.060
    % 10 3.708 6.552
    % 100 2.696 3.984
    % inf 2.605 3.781
    % chi^2 approximation, 95% fractile
    if strcmp(type,'chi') == 1
        semiaxes = sqrt(7.815*diag(dDxhat));
        title1 = '95% confidence ellipsoid, \chi^2 approx.';
    elseif strcmp(type,'F') == 1
        semiaxes = sqrt(3*6.591*diag(dDxhat));
        title1 = '95% confidence ellipsoid, F approx.';
    else
        error('drawConfidenceEllipsoid: wrong given type!');
        return;
    end
    h = figure;
    ellipsoidrot(0,0,0,semiaxes(1),semiaxes(2),semiaxes(3),vDxhat);
    axis equal
    xlabel(strcat('x [', metrics{1}, ']')); 
    ylabel(strcat('y [', metrics{2}, ']')); 
    zlabel(strcat('z [', metrics{3}, ']'));
    title(title1);
    view(37.5,15);
end