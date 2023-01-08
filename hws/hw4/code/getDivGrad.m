%% PROBLEM 2 - YOUR CODE HERE TO COMPUTE DIVERGENCE AND GRADIENT
% Div (|V| x 3|F| sparse) maps face-based vector fields to their vertex-based divergence
% Grad (3|F| x |V| sparse) maps a scalar function to its face-based gradient
function [Div, Grad] = getDivGrad(data)

% Div = sparse(data.nv, 3 * data.nf);
% Grad = sparse(3 * data.nf, data.nv);

%% discrete gradient
tripI = reshape(1:3*data.nf, [3 data.nf])';
tripI = repmat(tripI, [1 3]);

tripJ = repmat(data.triangles, [1 3]);
tripJ = tripJ(:, [1 4 7 2 5 8 3 6 9]);

edges = data.vertices(data.triangles(:, [3 1 2]), :) - data.vertices(data.triangles(:, [2 3 1]), :);
edges = reshape(edges, [size(data.triangles) 3]);

normals = repmat(data.faceNormals, [1 3]);
normals = reshape(normals, [data.nf 3 3]);
normals = permute(normals, [1 3 2]);

Afinv = repmat(1./(2*data.triangleAreas), [1 9]);

tripV = cross(normals, edges, 3);
tripV = reshape(tripV, [data.nf 9]);
tripV = tripV(:, [1 4 7 2 5 8 3 6 9]);
tripV = tripV .* Afinv;

Grad = sparse(tripI, tripJ, tripV, 3*data.nf, data.nv);

%% integrated divergence
tripI = repmat(data.triangles, [1 3]); tripI = tripI(:, [1 4 7 2 5 8 3 6 9]);
tripJ = reshape(1:3*data.nf, [3 data.nf])'; tripJ = repmat(tripJ, [1 3]);

%% trigonometry
coss =-dot(edges(:, [2, 3, 1], :), edges(:, [3, 1, 2], :), 3);
sins = vecnorm(cross(edges(:, [2, 3, 1], :), edges(:, [3, 1, 2], :), 3), 2, 3);
cots = coss ./ sins;

e11 = reshape(edges(:, 3, :), [data.nf 3]); e12 =-reshape(edges(:, 2, :), [data.nf 3]);
e21 = reshape(edges(:, 1, :), [data.nf 3]); e22 =-reshape(edges(:, 3, :), [data.nf 3]);
e31 = reshape(edges(:, 2, :), [data.nf 3]); e32 =-reshape(edges(:, 1, :), [data.nf 3]);

tripV = zeros(data.nf, 9);
tripV(:, 1:3) = e11 .* cots(:, 3) + e12 .* cots(:, 2);
tripV(:, 4:6) = e21 .* cots(:, 1) + e22 .* cots(:, 3);
tripV(:, 7:9) = e31 .* cots(:, 2) + e32 .* cots(:, 1);

Div = 0.5*sparse(tripI, tripJ, tripV, data.nv, 3*data.nf);

end
%%% END HOMEWORK PROBLEM