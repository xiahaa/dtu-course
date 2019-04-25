   imColor=imread('./data/library2.jpg');
 figure;imshow(imColor);
 if size(imColor,3) == 3
     imGray=rgb2gray(imColor);
 else
     imGray = imColor;
 end
 im = im2double(imGray);

 % here starts your code
 ......
 ......
 ......

 [row,col] = nonmaxsuppts(C,'radius', 2);

 % plot result
 figure
 img=imshow(imColor),title('my-Harris'),
 hold on
 plot(col,row, 'ro','MarkerSize',10),
 hold off
