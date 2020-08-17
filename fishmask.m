
    function BW = fishmask(A)
 %convert image to 2D gray scale 
    %A2d = rgb2gray(A);
 %use edge function to obtain threshold value of image 
    [~,threshold] = edge(A,'canny');
    fudgeFactor = 0.5;
 %use edge to obtain mask with tuned threshold
    Aedge = edge(A,'canny',threshold*fudgeFactor);
 % create convex hull to join up holes in image
    Aconv = bwconvhull(Aedge,'objects');
 % fill gaps by dilating canny image, create two perpendicular linear structing elements
     se90 = strel('line',3,90);
     se0 = strel('line',3,0);
     Adil = imdilate(Aconv,[se90 se0]);
 % fill interior gaps
    Afill = imfill(Adil,'holes');
% isolate things touching the border aka the fish and remove everything
% else
    Anobord = imclearborder(Afill,4);
    Amidd = Afill - Anobord;
%getting rid of small objects with under 8000 pixels (noise)
    Ano_noise = bwareaopen(Amidd, 8000);
% smooth image by eroding 
    seD = strel('rectangle',[1 2]);
    Aeroded = imerode(Ano_noise,seD);
%blur the image and then threshold at 0.5
    width = 15;
    kernel = ones(width) / width^2;
    blurryImage = conv2(double(Aeroded), kernel, 'same');
    Ablurry = blurryImage > 0.5;
%fill holes again 
    BW = imfill(Ablurry, 'holes');
padH = 000; % vertical padding
padW = 200;% horizontal padding
BW = padarray(BW,[padH, padW]);
%show og image next to mask 
    %imshow(mask)
    %title('Image Mask');
    end 
