function vectorDiffusion()

addpath('utils/');

% load mesh
meshfile = 'meshes/moomoo.off';
[X,T] = readOff(meshfile);
data = MeshData(X,T);
CL = buildConnectionLaplacian(data);

% build vector valued signal: delta vector heat source at highest Y axis point.
dt = sqrt(mean(data.edgeLengths));
signal = zeros(3*data.nv,1);
[~,maxind] = max(data.vertices(:,2));
t1 = cross(data.vertNormals(maxind,:),randn(3,1)'); t1 = t1/norm(t1);
signal((maxind-1)*3 + (1:3)) = t1;

%% PROBLEM 3(b) - VECTOR HEAT METHOD
% todo: compute diffused signal for various dt and normalize
% phi = zeros(data.nv, 3);

%% scalar Laplacian
L =-data.cotLaplacian;
A = data.massMatrix;

%% vector mass matrix
AV = repmat(data.vertexWeights, [1 3])';
AV = reshape(AV, 1, []);
AV = sparse(1:3*data.nv, 1:3*data.nv, AV);

%% vector heat flow Y(t)
Yt = (AV-dt*CL) \ signal;
Yt = reshape(Yt, [3 data.nv])';
Yt = Yt ./ vecnorm(Yt, 2, 2);

%% scalar heat flow u(t)
u0 = zeros(data.nv, 1); u0(maxind) = norm(t1);
ut = (A-dt*L) \ u0;

%% scalar heat flow phi(t)
phi0 = zeros(data.nv, 1); phi0(maxind) = 1.;
phit = (A-dt*L) \ phi0;

phi = (ut ./ phit) .* Yt;

%%% END HOMEWORK PROBLEM

%% Visualization
figure; trisurf(T, X(:, 1), X(:, 2), X(:, 3), 'FaceColor', 'k', 'LineStyle', 'none'); hold on;
q2 = quiver3(X(:,1),X(:,2),X(:,3),phi(:,1),phi(:,2),phi(:,3), 'r', 'LineWidth', 1);
q1 = quiver3(X(maxind,1),X(maxind,2),X(maxind,3),t1(1),t1(2),t1(3), 'c', 'LineWidth', 2, 'AutoScaleFactor', 3);
view(3); axis image vis3d off;
title('Parallel Transported Vector Field');
legend([q1, q2], 'Source', 'Parallel transport');

end

%% PROBLEM 3(a)
%%% YOUR CODE HERE - implement the connection laplacian operator
% data: mesh data structure
% CL: |3V| x |3V| "sparse connection laplacian operator"
function CL = buildConnectionLaplacian(data)
%     CL = sparse(3 * data.nv, 3 * data.nv);
    tripI = repmat([(data.triangles-1)*3+1 (data.triangles-1)*3+2 (data.triangles-1)*3+3], [1 3]);
    tripI = reshape(tripI, [data.nf 9 3]);
    tripI = tripI(:, [1 4 7 2 5 8 3 6 9], :);
    tripI = repmat(tripI, [1 1 3]);
    
    tripJ = repmat([(data.triangles-1)*3+1 (data.triangles-1)*3+2 (data.triangles-1)*3+3], [1 3]);
    tripJ = reshape(tripJ, [data.nf 9 3]);
    tripJ = permute(tripJ, [1 3 2]);
    tripJ = tripJ(:, :, [1 4 7 2 5 8 3 6 9]);
    tripJ = repmat(tripJ, [1 3 1]);
    
    tripV = zeros(data.nf, 9, 9);

    %% edges
    edges = data.vertices(data.triangles(:, [3 1 2]), :) - data.vertices(data.triangles(:, [2 3 1]), :);
    edges = reshape(edges, [size(data.triangles) 3]);
    
    %% trigonometry
    coss =-dot(edges(:, [2, 3, 1], :), edges(:, [3, 1, 2], :), 3);
    sins = vecnorm(cross(edges(:, [2, 3, 1], :), edges(:, [3, 1, 2], :), 3), 2, 3);
    cots = coss ./ sins;
    
    for i=1:data.nf
        %% vertex index
        vi = data.triangles(i, 1);
        vj = data.triangles(i, 2);
        vk = data.triangles(i, 3);
        
        %% cotangent
        a = cots(i, 1);
        b = cots(i, 2);
        c = cots(i, 3);
    
        %% rotation
        ni = data.vertNormals(vi, :); ni = ni ./ norm(ni);
        nj = data.vertNormals(vj, :); nj = nj ./ norm(nj);
        nk = data.vertNormals(vk, :); nk = nk ./ norm(nk);
        
        rij = vrrotvec2mat(vrrotvec(nj, ni)); rji = rij';
        rjk = vrrotvec2mat(vrrotvec(nk, nj)); rkj = rjk';
        rki = vrrotvec2mat(vrrotvec(ni, nk)); rik = rki';
    
        tripV(i, 1:3, 1:3) = (b+c)*eye(3); tripV(i, 1:3, 4:6) = -c*rij;       tripV(i, 1:3, 7:9) = -b*rik;
        tripV(i, 4:6, 1:3) = -c*rji;       tripV(i, 4:6, 4:6) = (c+a)*eye(3); tripV(i, 4:6, 7:9) = -a*rjk;
        tripV(i, 7:9, 1:3) = -b*rki;       tripV(i, 7:9, 4:6) = -a*rkj;       tripV(i, 7:9, 7:9) = (a+b)*eye(3);
    end

    tripI = reshape(tripI, data.nf, []);
    tripJ = reshape(tripJ, data.nf, []);
    tripV = reshape(tripV, data.nf, []);
    
    CL = sparse(tripI, tripJ, tripV, 3*data.nv, 3*data.nv);
    CL = -0.5*CL;

end
%%% END HOMEWORK PROBLEM

