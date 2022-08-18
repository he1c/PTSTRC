function index_set = patch_group(point, radius, search_rad, search_gap, ...
                                            img, THRESHOLD, mode)
% mode: 1: distance < THRESHOLD; 2: the first THRESHOLD patches


[row,col] = size(img);

if point(1) < radius+1 || point(1) > row-radius || point(2) < radius+1 || point(2) > col-radius
    error('Out of range');
end

s_patch = point2roi(point, img, radius);

x(1) = max([radius+1,point(1)-search_rad]);
x(2) = min([row-radius,point(1)+search_rad]);
x_range = x(1):search_gap:x(2); x_range_len = length(x_range);
y(1) = max([radius+1,point(2)-search_rad]);
y(2) = min([col-radius,point(2)+search_rad]);
y_range = y(1):search_gap:y(2); y_range_len = length(y_range);

point_set = zeros( x_range_len * y_range_len , 2);
for k = 1:y_range_len
    point_set( (k-1)*x_range_len+1 : k*x_range_len, :) = [x_range; linspace(y_range(k),y_range(k),x_range_len)]';
end

if mode == 1
    point_set_valid = false(size(point_set,1),1);
    for i = 1:size(point_set,1)
        p_patch = point2roi(point_set(i,:), img, radius);
        dif = sum(sum((s_patch-p_patch).^2));
        disp(dif);
        if dif < THRESHOLD
            point_set_valid(i) = true;
        end
    end

    index_set = point_set(point_set_valid,:);

else
    point_set_dif = zeros(size(point_set,1),1);
    for i = 1:size(point_set,1)
        p_patch = point2roi(point_set(i,:), img, radius);
        point_set_dif(i) = sum(sum((s_patch-p_patch).^2));
    end
    [~,sort_index] = sort(point_set_dif,'ascend');
    if point_set(sort_index(1),:) == point
        index_set = point_set(sort_index(2:THRESHOLD+1),:);
    else
        index_set = point_set(sort_index(1:THRESHOLD),:);
    end
end

