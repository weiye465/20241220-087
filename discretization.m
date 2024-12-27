% 结构离散化
n_elements = 20;           % 单元数量
element_length = L/n_elements;  % 单元长度
nodes = n_elements + 1;    % 节点数量

% 建立节点坐标
z_coords = linspace(-h-6*D, 0, nodes);  % 从底部(-h-6D)到水面(0) 