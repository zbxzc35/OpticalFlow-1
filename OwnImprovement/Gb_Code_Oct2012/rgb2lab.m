function im = rgb2lab(im)

cform = makecform('srgb2lab');
im = applycform(im,cform);
end