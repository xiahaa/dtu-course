   clc;close all;clear all;

 % generate 500 points
   [X, lineTrue] = gen_line_data(500);

 % inhomogeneous to homogeneous
   ph = tohomo(X);

   % your code starts here, do RANSAC line fitting,
   % the result is the parameters of the fitted line, i.e. params(1)*x+params(2)*y+params(3)=0
   params = xxxx(ph);

   % Results Visualization
   figure;
   hold on
   xmin = min(X(1,:));
   xmax = max(X(1,:));
   xx = linspace(xmin,xmax,100);
   yy1 = -(lineTrue(1).*xx+lineTrue(3))./(lineTrue(2)+1e-6);
   yy2 = -(params(1).*xx+params(3))./(params(2)+1e-6);
   plot(xx, yy1, 'k-.','LineWidth',1);
   plot(xx, yy2, 'm--','LineWidth',2);
   xlabel('x')
   ylabel('y')
   title('RANSAC results for 2D line estimation')
   axis equal tight
   set(gca,'FontName','Arial','FontSize',20);
