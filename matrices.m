% 组装总体矩阵
M = zeros(2*nodes, 2*nodes);  % 质量矩阵
K = zeros(2*nodes, 2*nodes);  % 刚度矩阵

% 单元矩阵
me = element_length/420 * [... 
    156, 22*element_length, 54, -13*element_length;
    22*element_length, 4*element_length^2, 13*element_length, -3*element_length^2;
    54, 13*element_length, 156, -22*element_length;
    -13*element_length, -3*element_length^2, -22*element_length, 4*element_length^2];

ke = EI/(element_length^3) * [...
    12, 6*element_length, -12, 6*element_length;
    6*element_length, 4*element_length^2, -6*element_length, 2*element_length^2;
    -12, -6*element_length, 12, -6*element_length;
    6*element_length, 2*element_length^2, -6*element_length, 4*element_length^2]; 