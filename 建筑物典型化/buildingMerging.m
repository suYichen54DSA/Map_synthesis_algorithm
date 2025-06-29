clc;
clear;
tic;
%% 算法输入 
N1 = 2000;
N2 = 60;
k = 1; %判定该超像素是否归属为T2区域的参数，超过k个sigma则判定该超参数为非T2区域，该参数可设置为小数
choose_centerORarea = 1; %选择以面积还是以中心为划分建筑物的基准
chooose_maxORmean = 1; %选择以一个超像素中最大面积建筑还是以建筑均值来作为典型化后建筑框的大小
building_data = imread("D:\Map synthesis algorithm\建筑物典型化\典型化测试数据\merg0303.png");
med_length = 10;
med_width = 10; %输入滤波核的长宽
items1 = 10; %第一次分割迭代的次数
items2 = 20; %第二次分割迭代的次数
%% 对带有建筑物的影像进行超像素分割
[L,numLabels] = superpixels(building_data,N1,"NumIterations",items1);
buildingNeighborSize = 10;
BW = boundarymask(L);
superpixel_types = 4 * ones(numLabels, 1);

%% 
se = strel('square', 5); % 设置形态学算子
building = building_data(:,:,1) == 0;    %生成建筑物掩膜
[label_building, num_labels] = bwlabel(building);
[buildingNeighbors] = buildingNeighbor(label_building,L,num_labels,buildingNeighborSize); %获得建筑物对象与其余超像素的关系

building_area = building.*L;   
building_label = unique(building_area); 
superpixel_types(building_label(2:end)) = 1;%用1表示建筑(B)
M_uni = unique (buildingNeighbors);
superpixel_types(M_uni(2:end)) = 2;%用2表示M
tabulated_B = tabulate(buildingNeighbors(buildingNeighbors ~= 0));
N_uni = tabulated_B(tabulated_B(:, 2) >= 2, 1);
superpixel_types(N_uni(2:end)) = 3;%用3表示N


%% 对BMN四类超像素进行边缘中值滤波平滑以生成合并区域T
building_BMN_all = zeros(size(building),'like',building); %首先需要对BMN的边界进行扩展，因此先生成一个逻辑矩阵
%创建一个逻辑数组，表示superpixel_types中是否等于'B'或'C'或'A'或'D'
is_BMN = find(ismember(superpixel_types, [1,2,3]));
%使用逻辑数组来获取所有'B'或'C'或'A'或'D'对应的point
index_BMN_all = ismember(L,is_BMN);
building_BMN_all(index_BMN_all) = 1;
mask_BMN = imfill(building_BMN_all,'holes'); %BMN区域内所有超像素的二值化掩膜
T = medfilt2(mask_BMN,[med_length med_width]);      %对其进行medfilt2操作，最终获得我们所需要的T区域
T_bw = bwmorph(T,'remove');         %获得其边界二值化图像
%%
%将红色边缘线绘制在原图像上，再进行超像素分割，而后再保留边缘线内的对分割结果，以此完成二次分割

building_L2 = ones(size(building),'like',building); 
building_L2 = building_L2.* ~T;


[L2,numLabels2] = superpixels(building_L2,N2,"NumIterations",items2);
mask_T2 = L2.* T;
Centroid = regionprops(L2,'Centroid');
% 舍弃与T区域有重合，但重合面积较小的区域
[uni_T2,~,label_index] = unique(mask_T2);
num_uni = histc(label_index, 1:numel(uni_T2));
num_uni = num_uni(2:end); %舍弃值为0的统计项
mean_num = mean(num_uni(2:end));
std_num = std(num_uni(2:end));
away_index = find(num_uni<k*mean_num-std_num); %根据正态分布的性质，舍弃小于2个sigma的值
away_index = [away_index+1;1];
uni_T2(away_index) = []; %将不符合的值舍弃
Centroid_T2 = Centroid(uni_T2);



Centroid_building = regionprops(label_building,'Centroid');%计算每个建筑物的中心坐标
% imshow(label_building)

[uni_building, ~, label_bui] = unique(label_building(:));
building_area = histc(label_bui, 1:numel(uni_building));   %获取每个建筑物对应的面积
[buildingBelongsArea] = buildingBelong(Centroid_building,L2,num_labels,uni_T2,label_building,building_area,buildingNeighborSize,choose_centerORarea); %获取各建筑物归属的超像素标记
% [buildingBelongs2] = buildingBelong(Centroid_building,L2,num_labels,uni_T2,label_building,building_area,buildingNeighborSize,0); %获取各建筑物归属的超像素标记

if chooose_maxORmean == 1 %以某超像素内最大面积为典型化后的面积
    buildingBelongsArea = max(buildingBelongsArea, [], 2);   %以最大面积为基准进行出图
elseif chooose_maxORmean == 0 %以某超像素内平均面积为典型化后的面积
    % 获取非零列的逻辑索引
    buildingBelongsArea(buildingBelongsArea == 0) = NaN; % 将0替换为NaN
    % 对非零列以行为基准取均值
    buildingBelongsArea = mean(buildingBelongsArea,2, 'omitnan');
end
%%
% 创建一个空白的 uint8 图像，用于存储填充后的结果
filled_image = uint8(zeros(size(building_data)));

% 定义不同类型超像素的颜色映射
color_map = [
    uint8([150, 218, 241]);  % 类型1的颜色
    uint8([148, 148, 148]);  % 类型2的颜色
    uint8([1, 153, 67]);     % 类型3的颜色
    uint8([251, 254, 0]);    % 类型4的颜色
];

% 根据超像素类型创建颜色矩阵
color_matrix = color_map(superpixel_types, :);

% 填充图像
filled_image = reshape(color_matrix(L(:), :), size(L, 1), size(L, 2), 3);
filled_image(BW) = 0; % 给超像素边界设置颜色

% 显示填充后的图像
figure
subplot(1, 3, 1)

% 创建图例的颜色示例
for label = 1:4
    color = color_map(label, :); % 获取对应标签的颜色
    plot(0, 0, 's', 'MarkerSize', 20, 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'DisplayName', char('A' + label - 1));
    hold on;
end

% 显示图例
lgd = legend('show');
set(lgd, 'NumColumns', 2); % 设置图例为2列

% 设置图例背景透明
set(lgd, 'Color', 'none');

% 调整图例中颜色块的长宽比为2:1
legend_color_boxes = findobj(lgd, 'Type', 'patch');
for i = 1:length(legend_color_boxes)
    set(legend_color_boxes(i), 'Position', get(legend_color_boxes(i), 'Position') + [0, 0, 0, -0.5]);
end

% 关闭坐标轴
axis off;
imshow(filled_image);
title("超像素分类结果图")

subplot(1, 3, 2)
T_building = building_data;
T_building(T_bw) = 0;
imshow(T_building);
title("T边界图");


subplot(1, 3, 3)

building_reg = 255*ones(size(building_data),'like',building_data); 
building_reg(T_bw) = 0;
building_reg(boundarymask(L2)~=0) = 122;
imshow(building_reg);
hold on



for i = 1:numel(Centroid_T2)
    if buildingBelongsArea(i) ~= 0  %若该超像素上所有建筑物均未被判定位于该超像素上，则不参与计算
        building_mask = (L2 == uni_T2(i)).*building;
        building_eccentricity = regionprops(building_mask,'Eccentricity');
        building_Orientation = regionprops(building_mask,'Orientation');
    
        centroid = Centroid_T2(i).Centroid;
        x = centroid(1); % x坐标
        y = centroid(2); % y坐标
        scatter(x, y, 5, 'b', 'filled'); % 这里的'50'是点的大小，'r'是点的颜色
 
        % 绘制矩形边框,根据前面的矩形归属来确定边框面积大小
        a = sqrt(buildingBelongsArea(i)/building_eccentricity.Eccentricity);
        b = a*building_eccentricity.Eccentricity; % 假设矩形的宽度
        angle_rad = -deg2rad(building_Orientation.Orientation); % 角度转换为弧度
    
        X = [0, a, a, 0, 0]; % 初始矩形边框顶点的x坐标
        Y = [0, 0, b, b, 0]; % 初始矩形边框顶点的y坐标
    
        % 旋转矩形坐标
        X_rotated = (X - a/2) * cos(angle_rad) - (Y - b/2) * sin(angle_rad) + x;
        Y_rotated = (X - a/2) * sin(angle_rad) + (Y - b/2) * cos(angle_rad) + y;
    
        h = plot(X_rotated, Y_rotated, 'b'); % 绘制旋转后的矩形边框，'b'是边框颜色
        
    end
end

title("建筑物典型化图");




%%
t1 = toc;
display(strcat('计算时间：',num2str(t1),'秒'));







