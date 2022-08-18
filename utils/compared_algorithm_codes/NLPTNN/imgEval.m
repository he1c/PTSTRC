function [RSE,PSNR,SSIM] = imgEval(img1, img2)

img1 = uint8(img1);
img2 = uint8(img2);

band = size(img1,3);

RSE = perfscore(img1, img2);
PSNR = psnr(img1, img2);
SSIM = 0;
for k = 1:band
    SSIM = SSIM + ssim_index(img1(:,:,k),img2(:,:,k));
end
SSIM = SSIM/band;