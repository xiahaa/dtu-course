
function P = triangulationIterLS(x1,P1,x2,P2)
    % initial weight
    w1 = 1;
    w2 = 1;
    
    iter = 1;
    maxIter = 50;
    oldCost = -1e6;
    
    while iter < maxIter
        % triangulate
        P = triangulationLS(x1,P1,x2,P2,1/w1,1/w2);
        % update
        nw1 = P1(3,:)*[P;1];
        nw2 = P2(3,:)*[P;1];
        % cost
        newCost = abs(nw1-w1)+abs(nw2-w2);
        
        if abs(newCost-oldCost) < 1e-6
            break;
        end
        oldCost = newCost;
        w1 = nw1; w2 = nw2;
        disp(iter);
        iter = iter + 1;
    end
    % final refine
    P = triangulationLS(x1,P1,x2,P2,1/w1,1/w2);
end