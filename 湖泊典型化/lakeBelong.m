function [lakeBelongsArea] = lakeBelong(Centroid_lake,T_result,num_labels,uni_T2,label_lake,lakeNeighborSize,choose,lake_data)  %判断lake和超像素的相邻关系
lakeBelongsArea = zeros(numel(uni_T2)-1,lakeNeighborSize);
if (choose == 1) %以中心为判定归属的基准
    center_lake = [];
    for i = 1:num_labels
        Centroid = Centroid_lake(i).Centroid;
        center_lake =  [center_lake;round(Centroid(1)),round(Centroid(2))];  %获取中心点
    end

    index = sub2ind(size(T_result),center_lake(:,2),center_lake(:,1));
    belongs = T_result(index); %通过各湖泊中心点坐标获取其对应超像素编号，对应的是T_result里面的编号

    % lake_data(index) = 255;
    % % color = label2rgb(T_result)
    % imshow(lake_data)
    for i = 2:numel(uni_T2)
        index_T2 = find(belongs == uni_T2(i));   %查询该湖泊存在于哪个超像素中
        lakeBelongsArea(i-1,1:numel(index_T2)) = index_T2;
    end
elseif (choose == 0) %以面积为判定归属的基准
    for i = 1:num_labels
        mask = gpuArray(T_result.*(label_lake == i));

        % mask_color = label2rgb(mask);  
        % imshow(mask_color)

        uni_lake = unique(mask);   %算出当前建筑物位于哪个超像素上
        j_total = [];
        for j = 2:numel(uni_lake)  %舍去0，从第2个被记录的超像素开始计算其面积
            j_total = [j_total;sum(mask == uni_lake(j),"all")];
        end 
        [~,max_j] = max(j_total);            %该结果中计算出的max_j即为对于该建筑物，包含面积最大的超像素序号
        index_T2 = find(uni_T2 == uni_lake(max_j+1));
        index_belongs = find(lakeBelongsArea(index_T2-1,:) == 0, 1);  %查找当前超像素对应行的最小为0的列数
        % lakeBelongsArea(index_T2,index_belongs) = lake_area(i+1);
        lakeBelongsArea(index_T2-1,index_belongs) = i;

    end

end