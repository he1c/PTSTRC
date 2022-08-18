clc 
clear;
close all;

addpath(genpath('F:\tensor\algorithm\gtnn'));
addpath(genpath('F:\tensor\algorithm\LRTC'));
addpath(genpath('F:\tensor\algorithm\LG'));
addpath(genpath('F:\tensor\algorithm\bayes'));
addpath(genpath('F:\tensor\algorithm\LRTC_TV'));
addpath(genpath('F:\tensor\algorithm\NLR\NLR_CS'));

% filename='airplane';
% filename='baboon';
% filename='barbara';
% filename='facade';
% filename='house';
filename='lena';
% filename='peppers';
% filename='sailboat';

filepath = 'F:\tensor\data\';
savePath = 'F:\tensor\algorithm\Patch\re\';

ObsRatio = [0.2,0.4,0.5,0.6,0.8];
radius = 3;
search_rad = 30;
search_gap = 1;
THRESHOLD = 5;
para.rho = 0.01; para.maxItr = 300; para.alpha = [1,1,1]; para.maxRank = 100;

RSE = zeros(length(ObsRatio),8);
PSNR = zeros(length(ObsRatio),8);
SSIM = zeros(length(ObsRatio),8);

img_path = [filepath,filename,'.bmp'];
img = imread(img_path);

img = double(img);
dim = size(img);

for fn =1:length(ObsRatio)

    Observ = mask_mat(2, dim, ObsRatio(fn));
    mask = ~Observ;

    img_out = zeros([dim, 8]);
    for i = 1:3
        img_out(:,:,i,1) = gray_inpaint(img(:,:,i), mask(:,:,i), radius, search_rad, search_gap, @t_SVD_new, THRESHOLD);
        img_out(:,:,i,2) = gray_inpaint(img(:,:,i), mask(:,:,i), radius, search_rad, search_gap, @HaLRTC_inpaint, THRESHOLD);
        img_out(:,:,i,3) = gray_inpaint(img(:,:,i), mask(:,:,i), radius, search_rad, search_gap, @FBCP_inpaint, THRESHOLD);
%         img_out(:,:,i,8) = NLR_inpaint(img(:,:,i), mask(:,:,i), mask(:,:,i), para);
    end

    img_init = img; img_init(mask) = mean(img(~mask));
    img_out(:,:,:,4) = t_SVD_inpaint(img, img_init, mask, para);
    img_out(:,:,:,5) = HaLRTC_inpaint(img, img_init, mask, para);
    img_out(:,:,:,6) = FBCP_inpaint(img, img_init, mask, para);
%     img_out(:,:,:,7) = TV_inpaint(img, img_init, mask, para);

    % reinforce for FBCP
    for k = 1:6
        tem = img_out(:,:,:,k);
        tem(~mask) = img(~mask);
        img_out(:,:,:,k) = tem;
    end
    
%     imwrite(img_out(:,:,:,1)/255, [savePath,filename,'_patch_tsvd_',num2str(ObsRatio(fn)),'.png'])
%     imwrite(img_out(:,:,:,2)/255, [savePath,filename,'_patch_halrtc_',num2str(ObsRatio(fn)),'.png'])
    imwrite(img_out(:,:,:,3)/255, [savePath,filename,'_patch_fbcp_',num2str(ObsRatio(fn)),'.png'])
%     imwrite(img_out(:,:,:,4)/255, [savePath,filename,'_raw_tsvd_',num2str(ObsRatio(fn)),'.png'])
%     imwrite(img_out(:,:,:,5)/255, [savePath,filename,'_raw_halrtc_',num2str(ObsRatio(fn)),'.png'])
    imwrite(img_out(:,:,:,6)/255, [savePath,filename,'_raw_fbcp_',num2str(ObsRatio(fn)),'.png'])
%     imwrite(img_out(:,:,:,7)/255, [savePath,filename,'_TV_',num2str(ObsRatio(fn)),'.png'])
%     imwrite(img_out(:,:,:,8)/255, [savePath,filename,'_NLR_',num2str(ObsRatio(fn)),'.png'])

    for k=1:8
        [RSE(fn,k),PSNR(fn,k),SSIM(fn,k)] = imgEval(img,img_out(:,:,:,k));
    end
end

RSE = RSE';
PSNR = PSNR';
SSIM = SSIM';