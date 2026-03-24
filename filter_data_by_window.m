function sub_data = filter_data_by_window(full_data, window_idx, total_points)
% filter_data_by_window: 按滚动窗口截取时间序列字段。

    if nargin < 3
        total_points = 96;
    end

    fields = fieldnames(full_data);
    sub_data = struct();
    for i = 1:length(fields)
        field = fields{i};
        if ismatrix(full_data.(field)) && size(full_data.(field), 2) == total_points
            sub_data.(field) = full_data.(field)(:, window_idx);
        elseif isvector(full_data.(field)) && length(full_data.(field)) == total_points
            sub_data.(field) = full_data.(field)(window_idx);
        else
            % 非时间序列字段直接复制
            sub_data.(field) = full_data.(field);
        end
    end
end
