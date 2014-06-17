%read image and make it gray

I1 = imread('car1.jpg');
I2 = imread('car2.jpg');

%%I1 to I2 flow field

winSize = 5;
ITER_NO = 3;
PYRE_NO = 3;
A = im2double(rgb2gray(I1));
B = im2double(rgb2gray(I2));
A = imresize( A, size(A) - mod( size(A), 2^(PYRE_NO-1) ) );
B = imresize( B, size(B) - mod( size(A), 2^(PYRE_NO-1) ) );


% construct pyramid
G = fspecial('gaussian',[3 3],1);
halfWindow = (winSize-1)/2;

Apyre = cell(PYRE_NO,1);
Bpyre = cell(PYRE_NO,1);
Apyre{1} = conv2( A, G, 'same' );
Bpyre{1} = conv2( B, G, 'same' );

for k = 2:PYRE_NO
    Apyre{k} = impyramid( Apyre{k-1}, 'reduce' );
    Bpyre{k} = impyramid( Bpyre{k-1}, 'reduce' );
end


% estimate flow field
for p = PYRE_NO:-1:1
    fprintf('Pyramid level: %d\n',p)
    
    A_p = imReflect( Apyre{p}, halfWindow);
    B_p = Bpyre{p};
    
    if (isequal(p,PYRE_NO))
        u = zeros(size( Apyre{p} ));
        v = zeros(size( Apyre{p} ));
        flag_ = 0;
    end
    
    for k = 1:ITER_NO
        fprintf('Pyramid no: %d, Iteration no: %d\n',p,k);
        B = imWarp( u, v, Bpyre{p} );
        
        B_ref = imReflect(B, halfWindow);
        [Ix Iy] = gradient( B_ref );
        
        
        H  = Hmatrix( Ix, Iy, halfWindow, 0.001 );
        
        It = A_p - B_ref;
        
        [us vs] = LKstep(It, Ix, Iy, H, halfWindow);
                 
        us = us(halfWindow+1:size(us,1)-halfWindow, halfWindow+1:size(us,2)-halfWindow);
        vs = vs(halfWindow+1:size(vs,1)-halfWindow, halfWindow+1:size(vs,2)-halfWindow);   
       
        u = u + us;
        v = v + vs;
    end
    
    if p ~= 1 
        u = 2 * imresize(u,size(u)*2,'bilinear');
        v = 2 * imresize(v,size(v)*2,'bilinear');
    end
    
end
Ir = imWarp(u,v,B);

% output gif
volume = zeros(size(I1,1),size(I1,2),size(I1,3),2);
volume(:,:,:,1) = im2double(I1);
volume(:,:,:,2) = im2double(I2);
frame2gif(volume,'input1.gif');
volume(:,:,:,2) = imWarp(u,v,im2double(I2));
frame2gif(volume,'input2.gif');

flow = zeros(size(u,1),size(u,2),2);
flow(:,:,1) = u;
flow(:,:,2) = v;
imflow = flowToColor(flow);
figure;imshow(imflow);
