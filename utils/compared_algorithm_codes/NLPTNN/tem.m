clc 
clear;
close all;

addpath(genpath('F:\tensor'));

filepath = 'F:\tensor\data\lena.bmp';
radius = 3;
search_rad = 30;
search_gap = 3;
THRESHOLD = 10;

img = imread(filepath);
img = double(img);
dim = size(img);

Observ = mask_mat(1, dim);
mask = Observ;

img_out = imgPrepro(img(:,:,1),mask(:,:,1));
img_out = img_out/255;
index_set = patch_group([15,20], radius, search_rad, search_gap, ...
                                            img_out, THRESHOLD, 2);
img_mis = img_out; img_mis(mask(:,:,1)) = 1;
for i=1:THRESHOLD
    subplot(2,5,i);
%     imshow(index_set2i,roi(index_set(i,:),img_out,radius))
    img_roi = img_mis(index_set(i,1)-radius:index_set(i,1)+radius,index_set(i,2)-radius:index_set(i,2)+radius,:);
    imshow(img_roi)
end

