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
filename='facade';
% filename='house';
% filename='lena';
% filename='peppers';
% filename='sailboat';

filepath = 'F:\tensor\data\';
savePath = 'F:\tensor\algorithm\Patch\re\';

ObsRatio = 0.3;
radius = 3;
search_rad = 30;
search_gap = 1;
THRESHOLD = 5;

maxItr = 5:5:500;

RSE = zeros(length(maxItr),3);
PSNR = zeros(length(maxItr),3);
SSIM = zeros(length(maxItr),3);

img_path = [filepath,filename,'.bmp'];
img = imread(img_path);

% img = double(img);
img = double(rgb2gray(img));

dim = size(img);

Observ = mask_mat(2, dim, 0.5);
mask = ~Observ(:,:,1);
img_pre = imgPrepro(img, mask);

for fn = 1:length(maxItr)
    para.maxItr = maxItr(fn);
    img_out = zeros([dim,3]);
    
    para.rho = 1;
    img_out(:,:,1) = t_SVD_new(img,img_pre,mask,para);
    img_out(:,:,2) = t_SVD_inpaint(img,img_pre,mask,para);
%     para.rho = 1;
%     img_out(:,:,3) = t_SVD_inpaint_2(img,img_pre,mask,para);
    
    
    for k=1:3
        [RSE(fn,k),PSNR(fn,k),SSIM(fn,k)] = imgEval(img,img_out(:,:,k));
    end
end

h1=plot(maxItr,PSNR(:,1),'-','LineWidth',2);hold on
h2=plot(maxItr,PSNR(:,2),'--','LineWidth',2);
str1='ours';str2='with $\mathcal{Q}^0=0$';
h = legend([h1,h2],str1,str2,'FontSize',50,'Location','southeast');
set(h,'Interpreter','latex')
xlabel('Iterations');ylabel('PSNR(dB)');