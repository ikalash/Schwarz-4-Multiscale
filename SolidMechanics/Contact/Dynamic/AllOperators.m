function [ P, R ] = AllOperators(X)

[num_levels, num_nodes_bottom] = size(X);

P = cell(num_levels, 1);
R = cell(num_levels, 1);

num_elements_bottom = num_nodes_bottom - 1;

num_elements_top = ...
    num_elements_bottom / 2^(num_levels - 1);

num_nodes_top = num_elements_top + 1;

num_nodes_level = num_nodes_top;

for level = 1 : num_levels - 1

    X_level = X(level, 1 : num_nodes_level)';
    [P_level, R_level] = Operators(X_level);

    P{level} = P_level;
    R{level} = R_level;

    num_nodes_level = 2 * (num_nodes_level - 1) + 1;

end
