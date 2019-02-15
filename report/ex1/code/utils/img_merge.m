function imc = img_merge(I1,I2)
% merge two image into one large image.
% Author: xiahaa@space.dtu.dk
    if size(I1)~=size(I2)
        I2 = imresize(I2,size(I1));
    end
    imc = cat(2,I1,I2);
end