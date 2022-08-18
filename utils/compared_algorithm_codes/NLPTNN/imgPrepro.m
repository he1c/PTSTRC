function img_out = imgPrepro(img, mask)
% mask : is_missing 
[m_col,m_row] = find(mask);

o_value = img(~mask);
[o_col,o_row] = find(~mask);

o_col = [0;size(img,2)+1;0;size(img,2)+1;o_col];
o_row = [0;0;size(img,1)+1;size(img,1)+1;o_row];
mean_value = repmat(mean(mean(img)),4,1);
o_value = [mean_value;o_value];

m_value = griddata(o_col,o_row,o_value,m_col,m_row);

img_out = img;
for i = 1:length(m_col)
    img_out(m_col(i),m_row(i)) = m_value(i);
    if isnan(m_value(i))
        img_out(m_col(i),m_row(i)) = max(img_out(m_col(i),m_row(i-1)),img_out(m_col(i-1),m_row(i)));
    end
end

end