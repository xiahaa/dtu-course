function X=Box3D()

nPerSide=100;	%The number of ponits per border in the box.
Xoff=-0.3; 		% the ofset allong the x-axis
Yoff=-0.3; 		% the ofset allong the y-axis
Zoff=3; 			% the ofset allong the z-axis

Range=(1:(nPerSide))/(nPerSide+1);
Zero=zeros(1,nPerSide);
One=ones(1,nPerSide);
X=[1 0 1 0 1 0 1 0;
   0 0 1 1 0 0 1 1;
   0 0 0 0 1 1 1 1
];
X=[X,[Range;Zero;Zero],[Zero;Range;Zero],[Zero;Zero;Range]];
X=[X,[Range;One;One],[One;Range;One],[One;One;Range]];
X=[X,[Range;Zero;One],[Zero;Range;One],[Zero;One;Range]];
X=[X,[Range;One;Zero],[One;Range;Zero],[One;Zero;Range]];
X=[X,[Range;One/3;Zero]];


X(1,:)=X(1,:)+Xoff;
X(2,:)=X(2,:)+Yoff;
X(3,:)=X(3,:)+Zoff;

