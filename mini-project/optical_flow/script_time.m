clc;close all;clear all;

baseDir = '../../data/optical_flow_data/';

filename = {'Army','Backyard','Basketball','Dumptruck','Evergreen','Grove','Mequon', ...
            'Schefflera','Teddy','Urban','Wooden','Yosemite'};
        
ts = zeros(size(filename,2), 2);

for fileid = 1:size(filename,2) 
    dir = strcat(baseDir, filename{fileid});

    fid = fopen(strcat(dir,'/time.txt'),'r');
    i = 1;
    while ~feof(fid)
        line = fgets(fid);
        txt = split(line);
        numtxt = str2double(txt);
        id = isnan(numtxt);
        ts(fileid,i) = numtxt(~id);
        i = i + 1;
    end
    fclose(fid);
end
figure
cmap = lines(2);
plot(ts(:,1),'-o','Color',cmap(1,:),'LineWidth',2);hold on;grid on;
plot(ts(:,2),'-d','Color',cmap(2,:),'LineWidth',2);
legend({'LK','HS'},'Location','northwest');
xticks([1:size(filename,2)]);
xlim([1 size(filename,2)]);
xticklabels(filename);
xtickangle(45);
title('Rumtime');
ylabel('Time: (s)');
set(gca,'FontName','Arial','FontSize',20);