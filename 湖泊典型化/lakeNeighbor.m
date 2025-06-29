function [lakeNeighbors] = lakeNeighbor(label_lake,label_superpixel,num_labels,lakeNeighborSize)  %判断lake和超像素的相邻关系
%%
BW = gpuArray(boundarymask(label_superpixel));
boundaries_superpixel = BW .* label_superpixel; %带标签的超像素边界
    
lakeNeighbors = zeros(num_labels,lakeNeighborSize);
SE = strel('disk',3); %设置形态学算子
for label = 1:num_labels
    mask = (label_lake == label);
    boundary = imdilate(mask,SE);
    boundary(mask) = 0; %取出边界
    neighbor = boundary.*boundaries_superpixel;
    uni = unique(neighbor)';
    lakeNeighbors(label,1:size(uni,2)-1) = uni(2:end);
end
