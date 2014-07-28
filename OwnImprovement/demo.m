%
warning('off'); %#ok<WNOFF>

%read image
I1 = im2double(imread('car1.jpg'));
I2 = im2double(imread('car2.jpg'));

[h,w,~] = size(I1);
width = min(h,w);
%%I1 to I2 flow field

winSize = 5;
minWidth = 20;
ratio = 0.75;
PYRE_NO = log(minWidth/width)/log(ratio);

baseSigma = 1/ratio-1;
n=round(log(0.25)/log(ratio));
nSigma=baseSigma*n;

Apyre = cell(PYRE_NO,1);
Bpyre = cell(PYRE_NO,1);
Apyre{1} = I1;
Bpyre{1} = I2;
for i = 2:PYRE_NO
    if(i<=n+1)
       sigma = baseSigma*(i-1);
       G = fspecial('gaussian',round(sigma*3),sigma);
       A = imfilter( I1, G, 'replicate' );
       B = imfilter( I2, G, 'replicate' );
       Apyre{i} = imresize(A,ratio^(i-1));
       Bpyre{i} = imresize(B,ratio^(i-1));
    else
       G = fspecial('gaussian',round(nSigma*3),nSigma);
       A = imfilter( Apyre{i-1-n}, G, 'replicate' );
       B = imfilter( Bpyre{i-1-n}, G, 'replicate' );
       %A = conv2( Apyre{i-n}, G, 'same' );
       %B = conv2( Bpyre{i-n}, G, 'same' );
       [nh,nw,~] = size(A);
       nwidth = min(nh,nw);
       rate=(ratio^(i-1))*width/nwidth;
       Apyre{i} = imresize(A,rate);
       Bpyre{i} = imresize(B,rate);
    end
end

%No need
%PYRE_NO = 3;
%A = imresize( I1, size(I1) - mod( size(I1), (1/ratio)^(PYRE_NO-1) ) );
%B = imresize( I2, size(I2) - mod( size(I2), (1/ratio)^(PYRE_NO-1) ) );


% construct pyramid
ITER_NO = 3;
halfWindow = (winSize-1)/2;


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
        [Ix,Iy] = gradient( B_ref );
        
        
        H  = Hmatrix( Ix, Iy, halfWindow, 0.001 );
        
        It = A_p - B_ref;
        
        [us,vs] = LKstep(It, Ix, Iy, H, halfWindow);
                 
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
