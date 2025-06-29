function [buildingNeighbors] = buildingNeighbor(label_building,label_superpixel,num_labels,buildingNeighborSize)  %判断building和超像素的相邻关系
%%
BW = gpuArray(boundarymask(label_superpixel));
boundaries_superpixel = BW .* label_superpixel; %带标签的超像素边界
    
buildingNeighbors = zeros(num_labels,buildingNeighborSize);
SE = strel('disk',3); %设置形态学算子
for label = 1:num_labels
    mask = (label_building == label);
    boundary = imdilate(mask,SE);
    boundary(mask) = 0; %取出边界
    neighbor = boundary.*boundaries_superpixel;
    uni = unique(neighbor)';
    buildingNeighbors(label,1:size(uni,2)-1) = uni(2:end);
end