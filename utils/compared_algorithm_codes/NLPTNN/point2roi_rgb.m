function img_roi = point2roi_rgb(point, img, radius)

img_roi = img(max(point(1)-radius,1):point(1)+radius,max(point(2)-radius,1):point(2)+radius,:);