function test_on_p

A = [zeros(3,3) eye(3); ...
     zeros(3,3) zeros(3,3)];
 
B = [zeros(3,3);eye(3)];

Kp = [15;15;30];
Kd = [12;12;10];

K = [-blkdiag(15,15,30) -blkdiag(12,12,10)];
 
Ac = A+B*K;


cvx_begin sdp
    variable P(6,6) symmetric
    minimize 0
    subject to
        (Ac'*P+P*Ac)<0
        P > 0
        P(1,1)>1
cvx_end


end