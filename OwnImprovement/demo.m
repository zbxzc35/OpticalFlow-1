%
warning('off'); %#ok<WNOFF>

%read image
I1 = im2double(imread('car1.jpg'));
I2 = im2double(imread('car2.jpg'));


% construct pyramid
[h,w,~] = size(I1);
width = min(h,w);
minWidth = 50;
ratio = 0.5;
PYRE_NO = ceil(log(minWidth/width)/log(ratio));

baseSigma = 1/ratio-1;
n=round(log(0.25)/log(ratio));
nSigma=baseSigma*n;

Apyre = cell(PYRE_NO,1);
Bpyre = cell(PYRE_NO,1);
Apyre{1} = I1;
Bpyre{1} = I2;
Ratiopyre = cell(PYRE_NO,1);
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
       [nh,nw,~] = size(A);
       nwidth = min(nh,nw);
       rate=(ratio^(i-1))*width/nwidth;
       Apyre{i} = imresize(A,rate);
       Bpyre{i} = imresize(B,rate);
    end
end

for i = 1:PYRE_NO-1
   Ratiopyre{i} = min(size(Apyre{i},1),size(Apyre{i},2))/min(size(Apyre{i+1},1),size(Apyre{i+1},2)); 
end


% windows search size
winSize = 5;
ITER_NO = 3;
halfWindow = (winSize-1)/2;
accThreshold = 0.0005;
fp = fopen('LKLog.txt','wt');
% estimate flow field
for p = PYRE_NO:-1:1
    fprintf(fp,'Pyramid level: %d\n',p);
    A_p = imReflect( Apyre{p}, halfWindow);
    A_p_Desaturate = imDesaturate(A_p);
    B_p = Bpyre{p};
    
    if (isequal(p,PYRE_NO))
        u = zeros(size( Apyre{p},1),size(Apyre{p},2));
        v = zeros(size( Apyre{p},1),size(Apyre{p},2));
        
    else
        %at this place, we won't use bilinear interpolation
        %instead, we should use convoluation to find the similarity
        
        u = imresize(u,[size(Apyre{p},1),size(Apyre{p},2)],'bilinear')*Ratiopyre{p};
        v = imresize(v,[size(Apyre{p},1),size(Apyre{p},2)],'bilinear')*Ratiopyre{p};
        
        %correct u and v
        A_p_conv = imReflect( Apyre{p}, 2*halfWindow);
        B_p_conv = imReflect(imWarp( u, v, Bpyre{p} ), halfWindow);
        for i = 1+2*halfWindow:size(A_p_conv,1)-2*halfWindow
            for j = 1+2*halfWindow:size(A_p_conv,1)-2*halfWindow
                %perform convoluation and search the minimal
                template = B_p_conv(i-halfWindow:i+halfWindow,j-halfWindow:j+halfWindow,:);
                matching = A_p_conv(i-2*halfWindow:i+2*halfWindow,j-2*halfWindow:j+2*halfWindow,:);
                C = conv2(matching(:,:,1),template(:,:,1),'same')...
                    +conv2(matching(:,:,2),template(:,:,2),'same')...
                    +conv2(matching(:,:,3),template(:,:,3),'same');
                [mini,minj] = find(C==min(C(:)));
                u(i,j) = mini-2*halfWindow;
                v(i,j) = minj-2*halfWindow;
            end
        end
    end
    
    for k = 1:ITER_NO
        fprintf(fp,'Pyramid no: %d, Iteration no: %d\n',p,k);
        
        %search the minimal correlation for each u and v
        B = imWarp( u, v, Bpyre{p} );
        B_ref = imReflect(B, halfWindow);
        B_ref_Desaturate = imDesaturate(B_ref);
        [Ix,Iy] = gradient( B_ref_Desaturate );
        H  = Hmatrix( Ix, Iy, halfWindow, 0.001 );
        It = A_p_Desaturate - B_ref_Desaturate;
        [us,vs] = LKstep(It, Ix, Iy, H, halfWindow);
        us = us(halfWindow+1:size(us,1)-halfWindow, halfWindow+1:size(us,2)-halfWindow);
        vs = vs(halfWindow+1:size(vs,1)-halfWindow, halfWindow+1:size(vs,2)-halfWindow);   
        u = u + us;
        v = v + vs;
        
        if(sum(sum(us.^2+vs.^2))/size(us,1)*size(us,2)<accThreshold)
            fprintf(fp,'Pyramid no: %d, Iteration no: %d\n break',p,k);
            break;
        end 
    end
    figure;imshow(B);
    pause;
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
