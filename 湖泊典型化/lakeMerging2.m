clc;
clear;
% close all;

tic;
%% 算法输入 
N1 = 2500;
N2 = 30;
k = 3;%筛去二次超像素分割后极小区域的参数,建议N2越小该值越小,对于小于20的N2,该值建议设为1~2.
lake_data = imread("D:\Map synthesis algorithm\建筑物典型化\典型化测试数据\Data03.png");
neighborSize = 20;
ratio = 0.75;%LSC法参数
choose1 = 1; %若选用1，则采用LSC法进行超像素分类,反之选用0则采用S
choose2 = 1;
items1 = 10; %第一次分割迭代的次数
items2 = 30; %第二次分割迭代的次数
%% 对带有湖泊的影像进行超像素分隔
if  choose1 == 1 
    L=double(LSC_mex(lake_data,N1,ratio));
    numLabels = max(L(:));
elseif choose1 == 0 
    [L,numLabels] = superpixels(lake_data,N1,'NumIterations',items1);
end
BW = boundarymask(L);
% figure
% imshow(lake_data)
% L = gpuArray(L);
% BW = gpuArray(BW); 
%获取各超像素间的邻接关系，这段跑得比较慢，可不用
% [Neighbors_SP,boundaries,neighbors_edge] = superpixelNeighbor(L,numLabels,neighborSize);  %获得各超像素相邻关系,用numLabels*neighborSize大小的矩阵表示
superpixel_types = zeros(numLabels,1); % 创建矩阵储存分类结果
se = strel('disk',3); %设置形态学算子

%% 由于所有分类均是以湖泊为参照物，因此在对超像素进行分类前，需先对湖泊进行分类。
% 对于一般问题可通过NDWI进行计算以提取水体，而在本数据中直接根据数据提取出水体
lake = lake_data(:,:,1) ~= 255; 
lakeNeighborSize = 20;
% % 使用bwlabel函数执行标签分析,这一段需要改
[label_lake, num_labels] = bwlabel(lake);
lakeNeighbors = lakeNeighbor(label_lake,L,num_labels,lakeNeighborSize);  %获得lake和超像素的相邻关系

%% 
lake_a = lake.*L;   
lake_label = unique(lake_a); 
superpixel_types(lake_label(2:end)) = 1;%用1表示A
%%
B_uni = unique (lakeNeighbors);
superpixel_types(B_uni(2:end)) = 2;%用2表示B
tabulated_B = tabulate(lakeNeighbors(lakeNeighbors ~= 0));
C_uni = tabulated_B(tabulated_B(:, 2) >= 2, 1);
superpixel_types(C_uni(2:end)) = 3;%用3表示C



%% 进行第二轮遍历，在该层遍历中先找出F
lake_ABC_all = zeros(size(lake),'like',lake); %首先需要对ABC的边界进行扩展，因此先生成一个逻辑矩阵
%创建一个逻辑数组，表示superpixel_types中是否等于'B'或'C'或'A'
is_ABC = find(ismember(superpixel_types, [1,2,3]));
%使用逻辑数组来获取所有'B'或'C'或'A'对应的point
index_ABC_all = ismember(L(:), is_ABC);

lake_ABC_all(index_ABC_all) = 1;
mask_ABC = imfill(lake_ABC_all,'holes'); %制作该掩膜用于将ABC内部区域扣除
%由于ABCDE五类元素均在ABC形成的缓冲区域中,因此可认为mask_ABC外的所有元素均可标记为F
mask_F = ~mask_ABC;
F_area = mask_F.*L;
F_uni = unique(F_area);
superpixel_types(F_uni(2:end)) = 6;%用6表示F

%% 进行第三轮遍历，由于D类标签在本题算法中是由B类标签转变而来，因此此时未被标记的超像素均为E类标签，同时对E类标签进行边界扩充
superpixel_types(find(superpixel_types == 0)) = 5;%用5表示E

% 若计算了superpixel，可尝试如下方法
% E_index = find(superpixel_types == 5);
% E_uni = Neighbors_SP(E_index);

% 若单纯只是为出图，可考虑如下算法
lake_E_all = zeros(size(lake),'like',lake);
is_E = find(ismember(superpixel_types, 5));
index_E_all = ismember(L,is_E);
lake_E_all(index_E_all) = 1;
lake_E_all = imfill(lake_E_all,'holes');
mask_E = imdilate(lake_E_all,se);
mask_E(index_E_all) = 0;
E_uni = mask_E.*L;

E_uni = unique(E_uni); %得到所有E类超像素周边的超像素序号
% is_B = find(superpixel_types) == 2; %找出属于B的点
% is_D = intersect(is_B,E_uni);  %寻找是B同时在D超像素周边的序号，即为D
% superpixel_types(is_D) = 4;  %用4表示D，即将在E周边的原本为B的超像素转化为D类超像素
superpixel_types(E_uni(2:end)) = 4;  %用4表示D，即将在E周边的原本为B的超像素转化为D类超像素

%% 对ABCD四类超像素进行边缘中值滤波平滑以生成合并区域T
lake_ABCD_all = zeros(size(lake),'like',lake); %首先需要对ABCD的边界进行扩展，因此先生成一个逻辑矩阵
%创建一个逻辑数组，表示superpixel_types中是否等于'B'或'C'或'A'或'D'
is_ABCD = find(ismember(superpixel_types, [1,2,3,4]));
%使用逻辑数组来获取所有'B'或'C'或'A'或'D'对应的point
index_ABCD_all = ismember(L,is_ABCD);
lake_ABCD_all(index_ABCD_all) = 1;
mask_ABCD = imfill(lake_ABCD_all,'holes'); %ABCD区域内所有超像素的二值化掩膜
T = medfilt2(mask_ABCD,[7 7]);      %对其进行medfilt2操作，最终获得我们所需要的T区域
% imshow(T)
T_bw = bwmorph(T,'remove');         %获得其边界二值化图像
%% 对T区域进行二次超像素划分
[label_T,T_n] = bwlabel(T); % 使用8连接方式
T_uint8 = uint8(T);
T_rgb = cat(3, T_uint8 * 255, T_uint8 * 255, T_uint8 * 255);

% figure;
Centroid_all = [];
T_result = zeros(size(lake_data,1),size(lake_data,2),'double');
last_num = 0; %用于保证对标签进行合并后，各区域不会出现重复
if  choose1 == 1 
    T_result=double(LSC_mex(T_rgb,N2,ratio));
    % color = label2rgb(T_result);
    % imshow(color)
    T_result = T_result.*T; 
elseif choose1 == 0 
    for i =1:T_n   
        mask = label_T == i;   %单独取出来
        [rows, cols] = find(mask);
        x_min = min(rows(:));
        x_max = max(rows(:));
        y_min = min(cols(:));
        y_max = max(cols(:));
        mask_lake = T_rgb(x_min:x_max,y_min:y_max,:);
        [L2,numLabels2] = superpixels(mask_lake,N2,'NumIterations',items2);
       
        % color = label2rgb(label_lake)
        % imshow(color)
        T_result(x_min:x_max,y_min:y_max,:) = T_result(x_min:x_max,y_min:y_max,:)+mask(x_min:x_max,y_min:y_max).*(L2+last_num).*lake(x_min:x_max,y_min:y_max);
        last_num = last_num+numLabels2;
    end
    % 获取唯一值
    unique_values = unique(T_result);
    
    % 为每个唯一值分配连续的整数值，使得T_result值连贯
    mapped_values = zeros(size(T_result));
    for i = 1:length(unique_values)
        mapped_values(T_result == unique_values(i)) = i-1;
    end
    T_result = mapped_values;
end

%% 处理极小区域
[uni_result, ~, label_result] = unique(T_result(:));
Tresult_area = histc(label_result, 1:numel(uni_result));   %获取每个湖泊对应的面积
Tresult_area = Tresult_area(2:end);
% 计算Tresult_area的均值和标准差
meanArea = mean(Tresult_area);
stdDevArea = std(Tresult_area);
% 找到面积小于均值减去2个标准差的区域索引
smallAreaIndices = find(Tresult_area < meanArea - k * stdDevArea);
% 在T_result中将对应区域赋值为0
for i = 1:numel(smallAreaIndices)
    currentLabel = uni_result(smallAreaIndices(i));
    T_result(T_result == currentLabel) = 0;
end
% 获取唯一值
uniqueValues = unique(T_result);
% 为每个唯一值分配连续的整数值，使得T_result值连贯
mappedValues = zeros(size(T_result));
for i = 1:length(uniqueValues)
    mappedValues(T_result == uniqueValues(i)) = i-1;
end
T_result = mappedValues;

%%


[uni_lake, ~, label_lak] = unique(label_lake(:));
lake_area = histc(label_lak, 1:numel(uni_lake));   %获取每个湖泊对应的面积


% Centroid_lake对应标签是i
Centroid_lake = struct('Centroid', []);
for i = 1:length(uni_lake)-1
    target_label = uni_lake(i+1);
    T_target = label_lake==target_label;
    Centroid_super = regionprops(T_target,'Centroid');
    Centroid_lake(end+1) = Centroid_super(1);
end
Centroid_lake(1) = [];
% Centroid_lake = regionprops(label_lake,'Centroid');
%% 遍历完全，展示分类结果


% 创建一个空白的 uint8 图像，用于存储填充后的结果
filled_image = uint8(zeros(size(lake_data)));

% 定义不同类型超像素的颜色映射
color_map = [
    uint8([150, 218, 241]);  % 类型1的颜色
    uint8([148, 148, 148]);  % 类型2的颜色
    uint8([1, 153, 67]);     % 类型3的颜色
    uint8([251, 254, 0]);    % 类型4的颜色
    uint8([1, 0, 250]);      % 类型5的颜色
    uint8([255, 255, 255])   % 类型6的颜色
];

% 根据超像素类型创建颜色矩阵
color_matrix = color_map(superpixel_types, :);

% 填充图像
filled_image = reshape(color_matrix(L(:), :), size(L, 1), size(L, 2), 3);
filled_image(BW) = 0; % 给超像素边界设置颜色

% 显示填充后的图像
figure
subplot(1, 2, 1)

% 创建图例的颜色示例
for label = 1:6
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

subplot(1, 2, 2)
T_lake = lake_data;
T_lake(T_bw) = 0;
imshow(T_lake);
title("T边界图");

lake_reg = 255*ones(size(lake_data),'like',lake_data); 
lake_reg(T_bw) = 0;
lake_reg(boundarymask(T_result)~=0) = 122;

T_area = bwlabel(T); % 将T区域标签化，便于后续的合并判定

% 获取边界
while true
    [uni_T2,~,~] = unique(T_result);
    % 找到问题：这里的lakeBelongsArea索引对应的是uni_T2从第二个开始的在T_result的序号，值对应为label_lake的序号，中心点也对应label_lake的序号
    [lakeBelongsArea] = lakeBelong(Centroid_lake,T_result,num_labels,uni_T2,label_lake,lakeNeighborSize,choose2,lake_data);
    lakeBelongsmax = zeros(size(lakeBelongsArea,1),1);
    Centroid_super_total = zeros(length(uni_T2)-1,3);
    for i = 1:length(uni_T2)-1
        target_label = uni_T2(i+1);
        T_target = T_result==target_label;

        Centroid_super = regionprops(T_target,'Centroid');
        if length(Centroid_super)>1
           se = strel('disk', 3);  % 为防止细碎的点对超像素中心产生影像,在这里对其进行开运算
           T_target = imopen(T_target, se);
           Centroid_super = regionprops(T_target,'Centroid'); %再划分中心
        end
        Centroid_super_total(i,1) = Centroid_super(1).Centroid(1);
        Centroid_super_total(i,2) = Centroid_super(1).Centroid(2);
        Centroid_super_total(i,3) = uni_T2(i+1);
    end

    %选出面积每个超像素中面积最大的湖泊
    for i = 1:size(lakeBelongsArea,1)
        nonzero = find(lakeBelongsArea(i,:) ~= 0);  %判断是否为空
        if ~isempty(nonzero)       
            area = lake_area(lakeBelongsArea(i,nonzero)+1);  
            [~,maxIndex] = max(area);     
            lakeBelongsmax(i) = lakeBelongsArea(i,maxIndex);  %对应label_lake的序号
        end
    end



    max_lake = zeros(size(label_lake),'double');
       
    for i = 1:length(Centroid_super_total)
        if lakeBelongsmax(i) ~= 0
            %对应label_lake的序号
            max_lake(label_lake == lakeBelongsmax(i)) = lakeBelongsmax(i);   %算出各超像素最大面积湖泊，并存为连续的label
        end
    end
    B = bwboundaries(max_lake);      %获取每个最大label对应的边界

    
    %% 最终版本的原因 %%
    %  他now_centroid是将lakeBelongsmax值顺序排序后所得的，所以下面需要手动改now_centroid使得其顺序与lakeBelongsmax有对应关系
    B = cell(0,1);
    now_centroid = struct('Centroid', []);
    for i = 1:length(lakeBelongsmax)
        if lakeBelongsmax(i) ~= 0
            target_label = lakeBelongsmax(i);
            T_target = max_lake==target_label;
            B = [B;bwboundaries(T_target)];
            Centroid_super = regionprops(T_target,'Centroid');
            now_centroid(end+1) = Centroid_super(1);
        end
    end
    now_centroid(1) = [];

    %设定boundary_cell_array的length
    max_size = max(cellfun(@length, B));
    boundary_array = zeros(max_size,2*length(B),'double');
    boundary_array_reshape = [];
    %初始化B和now_centroid的索引
    index = 0;
    nan_i = []; 
    for i = 1:length(Centroid_super_total)
        if lakeBelongsmax(i) ~= 0  %若该超像素上所有湖泊均未被判定位于该超像素上，则不参与计算 
            index = index+1; %索引自增

            % 计算平移向量
            translationVector = [Centroid_super_total(i,2) - now_centroid(index).Centroid(2),Centroid_super_total(i,1) - now_centroid(index).Centroid(1)]; 
            %问题:now_centroid对应的索引应该是lake的标号

            % translationVector = [0,0];
            % 进行边界平移
            boundary = unique(round(B{index} + translationVector),'rows');
            boundary_array(1:length(boundary),index*2-1:index*2) = boundary;
            subboundary = sub2ind(size(T_area),boundary(:,1),boundary(:,2));
            %储存3列信息，前两列为其归属点，第三列为归属的超像素编号
            boundary_array_reshape = [boundary_array_reshape;[boundary,repmat(index,size(boundary,1),1),T_area(subboundary)]];  %合并所有点坐标，并且在第三列加上其索引，表征其指向的湖泊编号
        else
            nan_i = [nan_i;i];
        end
    end
    for i = 1:length(nan_i)
        index_nan = find(Centroid_super_total(:,3) == nan_i(i));
        Centroid_super_total(index_nan,:) = [];  % 删除指定索引处的结构体
    end
    belong_area = [];
    for i = 1:numel(now_centroid)
        belong_area = [belong_area;i,T_area(round(now_centroid(i).Centroid(1)),round(now_centroid(i).Centroid(2)))];%储存超像素序号及对应的T区域序号为一行
    end
    %% 检查是否重合
    merged_labels = [];
    for i = 1:numel(now_centroid)
        cer_array = boundary_array_reshape(boundary_array_reshape(:,3)==i,:);
        other_array = boundary_array_reshape(~(boundary_array_reshape(:,3)==i),:);
        [~, idx_cer, idx_other] = intersect(cer_array(:,1:2), other_array(:,1:2), 'rows');  %求出当前湖泊与其他湖泊有无交集
        intersection = [cer_array(idx_cer,3), other_array(idx_other,3)]; %返回具有交集的湖泊对应的序号
        intersection_belong = [belong_area(intersection(:,1),2),belong_area(intersection(:,2),2);]; %取出重合超像素对应的T标签序号
        diff_index = ~(intersection_belong(:,1) == intersection_belong(:,2));  %取出intersection_belong中T标签序号不同的两个超像素的索引
        intersection(diff_index,:) = []; %在intersection中扣除其名单，即不对其进行合并，再进行下一步判定
        if ~isempty(intersection)  %若为非空，即判定两者存在交集
             merged_labels = [merged_labels; intersection]; 
        end
    end
    if (isempty(merged_labels))
        break;
    end
    merged_labels = unique(sort(merged_labels, 2), 'rows');
    idx = find(merged_labels(:, 1) == merged_labels(:, 2));
    merged_labels(idx, :) = [];

    
    % 将每两行中的标签合并为同一个标签
    for i = 1:size(merged_labels, 1)
        label1 = merged_labels(i, 1);
        label2 = merged_labels(i, 2);
    
        % 找到标签为 label1 和 label2 的像素
        mask_label1 = (T_result == label1);
        mask_label2 = (T_result == label2);
    
        % 计算标签为 label1 和 label2 的像素面积
        area_label1 = sum(mask_label1(:));
        area_label2 = sum(mask_label2(:));
    
        % 将面积小的类别合并到面积大的类别中
        if area_label1 <= area_label2
            T_result(mask_label1) = label2;  % 将标签为 label1 的像素合并到标签为 label2 的类别中
        else
            T_result(mask_label2) = label1;  % 将标签为 label2 的像素合并到标签为 label1 的类别中
        end
    end
    
    % 对标签进行处理以防止重复
    unique_labels = unique(T_result);
    unique_labels(unique_labels == 0) = [];  % 排除 0 标签
    unique_labels = sort(unique_labels);  % 对标签进行排序
    for i = 1:length(unique_labels)
        T_result(T_result == unique_labels(i)) = i;  % 对标签进行重新映射为连续的整数值
    end
end 

figure
lake_reg2 = 255*ones(size(lake_data),'like',lake_data); 
% lake_reg2(T_bw) = 0;
BW_T = boundarymask(T_result);
lake_reg2(T_bw | BW_T) = 0;
index_one = sub2ind(size(lake_reg2),boundary_array_reshape(:, 1),boundary_array_reshape(:, 2),ones(size(boundary_array_reshape, 1), 1) * 2);
index_two = sub2ind(size(lake_reg2),boundary_array_reshape(:, 1),boundary_array_reshape(:, 2),ones(size(boundary_array_reshape, 1), 1) * 1);
index_three = sub2ind(size(lake_reg2),boundary_array_reshape(:, 1),boundary_array_reshape(:, 2),ones(size(boundary_array_reshape, 1), 1) * 3);
lake_reg2([index_one;index_two;index_three]) = 0;
imshow(lake_reg2)
hold on;
for i = 1:length(Centroid_super_total)
       x = round(Centroid_super_total(i,1));  % 四舍五入得到最近的整数坐标
       y = round(Centroid_super_total(i,2));
       index = sub2ind(size(lake_reg), y, x, ones(size(x)) * 2);  % 计算索引
       plot(x, y, '.','MarkerSize', 20);  % 在图像上用红色星号标记中心点
end
hold off;
title("典型化结果图");

% 绘制边界
%%
t1 = toc;
display(strcat('计算时间：',num2str(t1),'秒'));
