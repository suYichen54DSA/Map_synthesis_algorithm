function [Neighbors,boundaries_superpixel,neighbors_edge] = superpixelNeighbor(L,numLabels,neighborSize)
%%
BW = gpuArray(boundarymask(L));
boundaries_superpixel = BW .* L; %带标签的边界
Neighbors = zeros(numLabels,neighborSize);
neighbors_edge = zeros(numLabels,neighborSize);
SE = strel('disk',3); %设置形态学算子
for label = 1:numLabels
    mask = (boundaries_superpixel == label);
    mask = imfill(mask,'holes');
    boundary = imdilate(mask,SE);
    boundary(mask) = 0; %取出边界
    neighbor = boundary.*boundaries_superpixel;
    [uni,~,label_index] = unique(neighbor);
    neighbor_edge = histc(label_index, 1:numel(uni));
    Neighbors(label,1:size(uni,1)-1) = uni(2:end);
    neighbors_edge(label,1:size(neighbor_edge,1)-1) = neighbor_edge(2:end);
end

