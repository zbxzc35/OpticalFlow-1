function [ B ] = imWarp( flowHor, flowVer, im)
%This function warps B towards A
B = zeros(size(im));
for i = 1:size(im,3)
    [x,y] = meshgrid(1:size(im,2),1:size(im,1));
    C = interp2(im(:,:,i), x+flowHor, y+flowVer, 'cubic');
    C(isnan(C)) = im(isnan(C));
    B(:,:,i) = C;
end
end
