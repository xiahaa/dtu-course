function varargout = navSolver(prs, sat_pos,verbose)
%% GNSS navigation solver: calculate antenna's position from at least 4
%% pseudoranges.
    

end


function x0 = DLT(prs, sat_pos)
%% use DLT to estimate an initial value
    %% initially, omit receiver clock error
    n = size(prs,1);
    v = 1:n;
    C = nchoosek(v,2);
    
    sat_pos_r = sat_pos(:,1).^2+sat_pos(:,2).^2+sat_pos(:,3).^2;
    
    dprs = prs(C(:,1)).^2 - prs(C(:,2)).^2 - sat_pos_r(C(:,1),1) + sat_pos_r(C(:,2),1);
    A = 2.*[sat_pos(C(:,2),1)-sat_pos(C(:,1),1), ...
            sat_pos(C(:,2),2)-sat_pos(C(:,1),2), ...
            sat_pos(C(:,2),3)-sat_pos(C(:,1),3)];
    x0 = lssolver(A,dprs);
    x0 = [x0;0];
end

function x0 = SOCP(prs,sat_pos)
%% estimate initial value using SOCP and relevant relaxation
    
end