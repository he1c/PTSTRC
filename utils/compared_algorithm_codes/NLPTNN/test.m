clc
clear;
close all;

% addpath(genpath('F:\tensor\algorithm\gtnn'));
% addpath(genpath('F:\tensor\algorithm\LRTC'));
% addpath(genpath('F:\tensor\algorithm\LG'));
% addpath(genpath('F:\tensor\algorithm\bayes'));

filepath = 'F:\tensor\data\peppers.bmp';
radius = 3;
search_rad = 30;
search_gap = 1;
THRESHOLD = 5;
method_num = 2;

img = imread(filepath);
img = double(img);
dim = size(img);

Observ = mask_mat(2, dim, 0.8);
mask = ~Observ;

img_out = zeros([dim, method_num]); img_pre = zeros(dim);
RSE = zeros(method_num,1);

for i = 1:3
    img_pre(:,:,i) = imgPrepro(img(:,:,i), mask(:,:,i));
    img_out(:,:,i,1) = gray_inpaint(img(:,:,i), mask(:,:,i), radius, search_rad, search_gap, @t_SVD_new, THRESHOLD);
    img_out(:,:,i,2) = gray_inpaint(img(:,:,i), mask(:,:,i), radius, search_rad, search_gap, @t_SVD_inpaint, THRESHOLD);
    %img_out(:,:,i,3) = gray_inpaint(img(:,:,i), mask(:,:,i), radius, search_rad, search_gap, @HaLRTC_inpaint, THRESHOLD);
    %img_out(:,:,i,4) = gray_inpaint(img(:,:,i), mask(:,:,i), radius, search_rad, search_gap, @FBCP_inpaint, THRESHOLD);
end

%para.rho = 0.01; para.maxItr = 300; para.alpha = [1,1,1]; para.maxRank = 100;

%img_init = img; img_init(mask) = mean(img(~mask));
%img_out(:,:,:,5) = t_SVD_inpaint(img, img_init, mask, para);
%img_out(:,:,:,6) = HaLRTC_inpaint(img, img_init, mask, para);
%img_out(:,:,:,7) = FBCP_inpaint(img, img_init, mask, para);

for k = 1:method_num
    RSE(k) = perfscore(img_out(:,:,:,k),img);
    subplot(2,4,k)
    imshow(img_out(:,:,:,k)/255)
end

RSE(method_num+1) = perfscore(img_pre,img);
