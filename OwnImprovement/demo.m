warning('off'); %#ok<WNOFF>

%read image
I1 = im2double(imread('car1.jpg'));
I2 = im2double(imread('car2.jpg'));


%construct pyramid
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

for i = PYRE_NO:-1:2
    Ratiopyre{i} = min(size(Apyre{i-1},1),size(Apyre{i-1},2))/min(size(Apyre{i},1),size(Apyre{i},2));
end
Ratiopyre{1} = 1;

%obtain the region of using convolution
ImageDifference = cell(PYRE_NO,1);
for i = 1:PYRE_NO
    ImageDifference{i} = rgb2gray(imsubtract(Apyre{i},Bpyre{i}));
end

addpath('./Gb_Code_Oct2012');
ImageEdgeSeg = cell(PYRE_NO,1);
for i = 1:PYRE_NO
    ImageEdgeSeg{i} = Gb_CSG( Bpyre{i} );
end

ImageEdgeCanny = cell(PYRE_NO,1);
for i = 1:PYRE_NO
    ImageEdgeCanny{i} = edge(rgb2gray(Bpyre{i}),'canny',0.5);
end

ImageEdgeGB = cell(PYRE_NO,1);
[nRows, nCols, aux] = size(I2);
imDiag = norm([nRows, nCols]);
wC = round(0.016*imDiag);
alpha_AB = 1.9;
for i = 1:PYRE_NO
    ImageEdgeCanny{i} = GbC_lambda( Bpyre{i}, wC, alpha_AB);
end

%search candidate for each layer
CandiateMask = cell(PYRE_NO,1);
for i = 1:PYRE_NO
    seg = ImageEdgeSeg{i};
    dif = ImageDifference{i};
    value1 = seg(find(seg~=0));
    value2 = dif(find(dif~=0));
    mask = (seg>median(value1))&(dif>median(value2));
    CandiateMask{i} = mask;
end

% windows search size
winSize = 5;
ITER_NO = 3;
halfWindow = (winSize-1)/2;
accThreshold = 0.0005;


% estimate flow field
for p = PYRE_NO:-1:1
    fprintf('Pyramid level: %d\n',p);
    A_p = imReflect( Apyre{p}, halfWindow);
    A_p_Desaturate = imDesaturate(A_p);
    B_p = Bpyre{p};
    
    if (isequal(p,PYRE_NO))
        u = zeros(size( Apyre{p},1),size(Apyre{p},2));
        v = zeros(size( Apyre{p},1),size(Apyre{p},2));
        
    else
        u = imresize(u,[size(Apyre{p},1),size(Apyre{p},2)],'bilinear')*Ratiopyre{p};
        v = imresize(v,[size(Apyre{p},1),size(Apyre{p},2)],'bilinear')*Ratiopyre{p};
        
        %correct u and v
        A_p_conv = imReflect( Apyre{p}, 2*halfWindow);
        B_p_conv = imReflect(imWarp( u, v, Bpyre{p} ), halfWindow);
        mask = CandiateMask{p};
        for i = 1:size(Apyre{p},1)
            for j = 1:size(Apyre{p},2)
                %add conditions here
                if(isOnEdge(mask,i,j))
                    %perform convoluation and search the minimal
                    template = B_p_conv(i:i+halfWindow*2,j:j+halfWindow*2,:);
                    matching = A_p_conv(i:i+2*halfWindow*2,j:j+2*halfWindow*2,:);
                    C = conv2(matching(:,:,1),template(:,:,1),'full')...
                        +conv2(matching(:,:,2),template(:,:,2),'full')...
                        +conv2(matching(:,:,3),template(:,:,3),'full');
                    [maxi,maxj] = find(C==max(C(:)),1,'first');
                    u(i,j) = maxi-3*halfWindow;
                    v(i,j) = maxj-3*halfWindow;
                end
            end
        end
    end
    
    for k = 1:ITER_NO
        fprintf('Pyramid no: %d, Iteration no: %d\n',p,k);
        
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
            fprintf('Pyramid no: %d, Iteration no: %d\n break',p,k);
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
