function img_roi = point2roi(point, img, radius)

img_roi = img(point(1)-radius:point(1)+radius,point(2)-radius:point(2)+radius);