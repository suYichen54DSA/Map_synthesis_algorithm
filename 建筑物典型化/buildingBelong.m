function [buildingBelongsArea] = buildingBelong(Centroid_building,L2,num_labels,uni_T2,label_building,building_area,buildingNeighborSize,choose)  %判断building和超像素的相邻关系
buildingBelongsArea = zeros(numel(uni_T2),buildingNeighborSize);
if (choose == 1) %以中心为判定归属的基准
    center_building = [];
    for i = 1:num_labels
        Centroid = Centroid_building(i).Centroid;
        center_building =  [center_building;round(Centroid(1)),round(Centroid(2))];  %获取中心点
    end

    index = sub2ind(size(L2),center_building(:,2),center_building(:,1));
    belongs = L2(index); %通过各建筑物中心点坐标获取其对应超像素编号
    for i = 1:numel(uni_T2)
        index_T2 = find(belongs == uni_T2(i));
        buildingBelongsArea(i,1:numel(index_T2)) = building_area(index_T2+1);
    end
elseif (choose == 0) %以面积为判定归属的基准
    for i = 1:num_labels
        mask = gpuArray(L2.*(label_building == i));

        % mask_color = label2rgb(mask);  
        % imshow(mask_color)

        uni_building = unique(mask);   %算出当前建筑物位于哪个超像素上
        j_total = [];
        for j = 2:numel(uni_building)  %舍去0，从第2个被记录的超像素开始计算其面积
            j_total = [j_total;sum(mask == uni_building(j),"all")];
        end 
        [~,max_j] = max(j_total);            %该结果中计算出的max_j即为对于该建筑物，包含面积最大的超像素序号
        index_T2 = find(uni_T2 == uni_building(max_j+1));
        index_belongs = find(buildingBelongsArea(index_T2,:) == 0, 1);  %查找当前超像素对应行的最小为0的列数
        buildingBelongsArea(index_T2,index_belongs) = building_area(i+1);
    end

end