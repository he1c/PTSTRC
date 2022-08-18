function img_out = inpaint_multi(imgm, maskm, radius, search_rad, search_gap, inpaint_algo, THRESHOLD)
% mask --- is_missing_mat

img=imgm(:,:,2);
mask=maskm(:,:,2);

[row,col] = size(img);

x(1) = radius+1; x(2) = row-radius;
x_range = x(1):radius:x(2); x_range_len = length(x_range);
if x_range(x_range_len) ~= x(2)
    x_range = [x_range,x(2)];
    x_range_len = x_range_len + 1;
end
y(1) = radius+1; y(2) = col-radius;
y_range = y(1):radius:y(2); y_range_len = length(y_range);
if y_range(y_range_len) ~= y(2)
    y_range = [y_range,y(2)];
    y_range_len = y_range_len + 1;
end

all_points = zeros( x_range_len*y_range_len, 2);
for k = 1:y_range_len
    all_points( (k-1)*x_range_len+1 : k*x_range_len, :) = [x_range; linspace(y_range(k),y_range(k),x_range_len)]';
end

all_points_num = size(all_points,1);
out_buffer = zeros(2*radius+1, 2*radius+1, THRESHOLD+1, all_points_num);
point_buffer = zeros(THRESHOLD, 2, all_points_num);
for i=1:1:3
    img_pre{i} = imgPrepro(imgm(:,:,i), maskm(:,:,i));
end
frame_set_all=zeros(all_points_num,THRESHOLD+1);
parfor k = 1:all_points_num
    s_img = point2roi(all_points(k,:),img,radius);
    s_mask = point2roi(all_points(k,:),mask,radius);
    s_init = point2roi(all_points(k,:),img_pre{2},radius);
    [index_set,frame_set] = patch_group_multi(all_points(k,:), radius, search_rad, search_gap, img_pre, THRESHOLD, 2); % mode 2
    g_tensor = zeros(size(s_img,1),size(s_img,2),THRESHOLD+1);
    g_mask = g_tensor;  g_init = g_tensor;
    g_tensor(:,:,1) = s_img; g_mask(:,:,1) = s_mask; g_init(:,:,1) = s_init;
    for j = 1:THRESHOLD
        index_img = point2roi(index_set(j,:),imgm(:,:,frame_set(j)),radius);
        index_mask = point2roi(index_set(j,:),maskm(:,:,frame_set(j)),radius);
        index_init = point2roi(index_set(j,:),img_pre{frame_set(j)},radius);
        g_tensor(:,:,j+1) = index_img;
        g_mask(:,:,j+1) = index_mask;
        g_init(:,:,j+1) = index_init;
    end
    g_tensor_out = inpaint_algo(g_tensor, g_init, g_mask);
    out_buffer(:,:,:,k) = g_tensor_out(:,:,:);
    point_buffer(:,:,k) = index_set;
    frame_set_all(k,:)=[2;frame_set];
end

for i = 1:all_points_num
    img_out(all_points(i,1)-radius:all_points(i,1)+radius,all_points(i,2)-radius:all_points(i,2)+radius) = out_buffer(:,:,1,i);
end

for i = 1:all_points_num
    for j = 1:THRESHOLD
        if frame_set_all(i,j+1)==2
            point_tem = point_buffer(j,:,i);
            img_out(point_tem(1)-radius:point_tem(1)+radius,point_tem(2)-radius:point_tem(2)+radius) = ...
                        (point2roi(point_tem,img_out,radius) + out_buffer(:,:,j+1,i))/2;
        end
    end
end

end

