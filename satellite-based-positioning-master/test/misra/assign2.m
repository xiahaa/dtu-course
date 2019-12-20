function assign2
    %% 1
    ca1 = genCA(1, 1023, 3, 6);
    %% 2
    ca2 = genCA(1024,2046, 3, 6);
    
    figure;
    plot(ca1,'b');
    hold on;grid on;ylim([0,3]);
    plot(ca2,'r')
    
    %% 3
    ca3 = genCA(1, 1023, 5, 7);
    
    %% 4
    ca4 = genCA(1, 1023, 1, 9);
    
    %% 5
    ca1new = (ca1 == 0)*1 + (ca1 == 1)*(-1);
    ca3new = (ca3 == 0)*1 + (ca3 == 1)*(-1);
    ca4new = (ca4 == 0)*1 + (ca4 == 1)*(-1);
    
    %% 6
    n = -1022:1022;
    auto_corr = zeros(numel(n),1);
    for i = 1:numel(n)
        shift = n(i);
        id = [0:1022]+shift;
        id = mod(id,1023) + 1;
        cashift = ca1new(id);
        auto_corr(i) = ca1new'*cashift;
    end
    figure
    plot(auto_corr);
    
    ca1s = genCA(201, 1023+200, 3, 6);
    ca1snew = (ca1s == 0)*1 + (ca1s == 1)*(-1);
    n = -1022:1022;
    auto_corr = zeros(numel(n),1);
    for i = 1:numel(n)
        shift = n(i);
        id = [0:1022]+shift;
        id = mod(id,1023) + 1;
        cashift = ca1new(id);
        auto_corr(i) = ca1snew'*cashift;
    end
    figure
    plot(auto_corr);
    
    n = -1022:1022;
    auto_corr = zeros(numel(n),1);
    for i = 1:numel(n)
        shift = n(i);
        id = [0:1022]+shift;
        id = mod(id,1023) + 1;
        cashift = ca1new(id);
        auto_corr(i) = ca3new'*cashift;
    end
    figure
    plot(auto_corr);
    
    cals = ca1new([350:end, 1:(350-1)]);
    ca3s = ca3new([905:end, 1:(905-1)]);
    ca4s = ca4new([75:end, 1:(75-1)]);
    casum = ca1s + ca3s + ca4s;
    casumn = casum + 4*randn(1023,1);
    n = 0:1022;
    auto_corr = zeros(numel(n),1);
    for i = 1:numel(n)
        shift = n(i);
        id = [0:1022]+shift;
        id = mod(id,1023) + 1;
        cashift = ca1new(id);
        auto_corr(i) = casumn'*cashift;
    end
    figure
    plot(auto_corr);
end

function ca = genCA(epochs, epoche, s1, s2)
    G = ones(2,10);
    PNR = [];

    id = epochs:epoche;
    id = mod(id, 1023);
    id(id == 0) = 1;
    
    for i = 1:1023
        g1 = G(1,end);
        g2 = xor(G(2,s1)==1,G(2,s2)==1);
        g3 = xor(g1==1,g2==1);
        PNR = [PNR;g3];
        newg1 = xor(g1==1,G(1,3)==1);
        newg2 = xor(xor(xor(xor(xor(G(2,2)==1,G(2,3)==1),G(2,6)==1),G(2,8)==1),G(2,9)==1),G(2,10)==1);
        G(1,:) = [newg1,G(1,1:9)];
        G(2,:) = [newg2,G(2,1:9)];
    end
    ca = PNR(id);
end