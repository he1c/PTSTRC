clc
close all;

addpath(genpath('F:\tensor\algorithm\LG'));

% filepath = 'F:\tensor\data\lena.bmp';
% img = imread(filepath);
% img = double(rgb2gray(img))/255;

K=80;S=10;
img = zeros(K);
for i=-K+1:S:K-1
    for j = 0:1:S/2-1
        img = img+diag(linspace(0,0,K-abs(i+j)),i+j);
    end
    for j = S/2:1:S-1
        img = img+diag(linspace(0.8,0.8,K-abs(i+j)),i+j);
    end
end
% for i=1:6:43
%     img(i:i+1,:) = 0.4;
%     img(i+2:i+3,:) = 0.8;
%     img(i+4:i+5,:) = 0;
% end

radius = 5;
THRESHOLD = 4;
point = [25,25];
dim = size(img);
mask = ~mask_mat(2, dim, 0.8); mask = mask(:,:,1);
img_pre = imgPrepro(img, mask);

index_set = patch_group(point, radius, 10, 3, img_pre, THRESHOLD, 2);
index_set_ori = patch_group(point, radius, 10, 1, img, THRESHOLD, 2);

figure
for i=1:THRESHOLD
    subplot(3,THRESHOLD,i);
    imshow(point2roi(index_set(i,:),img_pre,radius))
    title(num2str(index_set(i,:)))
    subplot(3,THRESHOLD,i+THRESHOLD);
    imshow(point2roi(index_set(i,:),img,radius))
    subplot(3,THRESHOLD,i+2*THRESHOLD);
    imshow(point2roi(index_set_ori(i,:),img,radius))
    title(num2str(index_set_ori(i,:)))
end

imwrite(img_pre,'img_pre.bmp');
imwrite(img,'img.bmp');
img_mis = img;
img_mis(mask) = 1;
imwrite(img_mis,'img_mis.bmp');