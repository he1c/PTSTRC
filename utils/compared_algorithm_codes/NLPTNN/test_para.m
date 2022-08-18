clc
clear;
close all;

addpath(genpath('F:\tensor'));

filepath = 'F:\tensor\data\peppers.bmp';
radius = [6,8,10];
search_rad = 30;
search_gap = 3;
THRESHOLD = [5,10,15];

img = imread(filepath);
img = double(img)/255;
dim = size(img);

Observ = mask_mat(2, dim, 0.3);
mask = ~Observ;

RSE = zeros(3,3);
for k1 = 1:3
    for k2 = 1:3
        img_out = zeros(dim);
        tic
        for i = 1:3
            img_out(:,:,i) = gray_inpaint(img(:,:,i), mask(:,:,i), radius(k1), search_rad, search_gap, @t_SVD_inpaint,  THRESHOLD(k2));
        end
        toc
        RSE(k1,k2) = perfscore(img_out,img);
        subplot(1,2,1)
        imshow(img_out)
        subplot(1,2,2)
        img_out(mask) = 1;
        imshow(img_out)
    end
end


