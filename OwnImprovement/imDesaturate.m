function im = imDesaturate(im)


if(size(im,3)==3)
    im = im(:,:,1)*0.299+im(:,:,2)*0.587+im(:,:,3)*0.114;
end

end