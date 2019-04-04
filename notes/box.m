clear
close all
P=Box3D;
plot3(P(1,:),P(2,:),P(3,:),'.'),
axis equal
axis([-1 1 -1 1 -1 5])
xlabel('x')
ylabel('y')
zlabel('z')
