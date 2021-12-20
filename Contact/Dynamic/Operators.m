function [ P, R, H, K, M ] = Operators(X)

num_nodes_level = length(X);
num_elements_level = num_nodes_level - 1;
num_nodes_below = 2 * num_elements_level + 1;

% For this simple case, we can pre-compute the restriction
% and prolongation operators.
H = zeros(num_nodes_level, num_nodes_level);
K = zeros(num_nodes_level, num_nodes_below);
M = zeros(num_nodes_below, num_nodes_below);

for i = 1 : num_elements_level

  l = X(i + 1) - X(i);
  H_level = l * [1/3, 1/6; 1/6, 1/3];
  K_level = l * [5/24, 1/4, 1/24; 1/24, 1/4, 5/24];
  M_level = l * [1/6, 1/12, 0; 1/12, 1/3, 1/12; 0, 1/12, 1/6];

  H(i : i + 1, i : i + 1) = ...
      H(i : i + 1, i : i + 1) + H_level;

  K(i : i + 1, 2*i - 1 : 2*i + 1) = ...
      K(i : i + 1, 2*i - 1 : 2*i + 1) + K_level;

  M(2*i - 1 : 2*i + 1, 2*i - 1 : 2*i + 1) = ...
      M(2*i - 1 : 2*i + 1, 2*i - 1 : 2*i + 1) + M_level;

end

P = inv(M) * K';
R = inv(H) * K;
