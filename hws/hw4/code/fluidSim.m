%% PROBLEM 4
%%% a very basic intrinsic surface fluid simulation 
% function fluidSim()

%% clear
clear; clc;

addpath('utils/');

% load mesh. 
meshfile = 'meshes/spherers.off';
[X,T] = readOff(meshfile);
data = MeshData(X,T);
[Div, Grad] = getDivGrad(data);
meanEdgeLength = mean(data.edgeLengths);
dt = 0.0005;

% initialize velocity field and timestep
V = [zeros(data.nf, 1), -data.faceCenters(:, 3), data.faceCenters(:, 2)];
V = V - sum(V .* data.faceNormals, 2) .* data.faceNormals;
V = V ./ vecnorm(V, 2, 2);
V = (10 * meanEdgeLength) .* V;

% initialize color field
color = hsv2rgb([cos(3 * data.vertices(:, 3)).^2, ones(data.nv, 2)]);% + sin(10 * data.vertices(:, 2));

% initialize source region
sourceInd = find(max(data.faceCenters(:, 3)) - data.faceCenters(:, 3) < 0.05);

% visualize starting conditions
figure; hold all; rotate3d on; axis image vis3d off; set(gcf,'color','w');
ptc = patch('vertices', data.vertices,'faces',data.triangles,'edgecolor','none','FaceVertexCData', color);
shading interp;
qvr = quiver3(data.faceCenters(:,1),data.faceCenters(:,2),data.faceCenters(:,3),V(:,1),V(:,2),V(:,3),'b','linewidth',2);
scatter3(data.faceCenters(sourceInd,1),data.faceCenters(sourceInd,2),data.faceCenters(sourceInd,3),'b');

decL = decomposition(data.cotLaplacian);

% simulation loop
maxiters = 100000;
t = 0;
for i=1:maxiters
    sourceForce = (0.2 * meanEdgeLength / dt) * [cos(0.2 * pi * t) sin(0.2 * pi * t) 0];
    sourceForce = sourceForce - sum(sourceForce .* data.faceNormals, 2) .* data.faceNormals;

    %%% PROBLEM 4(d) - YOUR CODE HERE
    % integrate source force
    V = V + sourceForce * dt;
    
    %%% PROBLEM 4(e) - YOUR CODE HERE
    % implement pressure projection. Make sure V is divergence free
    VV = reshape(V', [3*data.nf, 1]);
    zeta =-decL \ (Div * VV);

    V = V - reshape(Grad * zeta, [3, data.nf])';
    
    % visualize velocity field
    if mod(i, 66) == 0
        qvr.UData = V(:, 1);
        qvr.VData = V(:, 2);
        qvr.WData = V(:, 3);
        ptc.FaceVertexCData = color;
        drawnow;
        pause(0.01);
    end
    
    %%% PROBLEM 4(a) - YOUR CODE HERE
    % construct Dv operator
    Dv = getDirectionalDerivative(data, V, Grad);
    
    %%% PROBLEM 4(b) - YOUR CODE HERE
    % advect color field
    TV = expMap(Dv, -dt);
    color = TV * color;

    %%% PROBLEM 4(c) - YOUR CODE HERE
    % implement fluid self-advection
    
    %%% END HOMEWORK PROBLEM
    
    t = t + dt;
end

% end

function Dv = getDirectionalDerivative(data, Vt, Grad)

%% directional derivative
tripI = repmat((1:data.nf)', 1, 3);
tripJ = repmat((((1:data.nf) - 1)*3)', 1, 3) + (1:3);

Dv = sparse(tripI, tripJ, Vt, data.nf, 3*data.nf) * Grad;

%% face to vertex
Dv = data.FtoV * Dv;
end

function T = expMap(D, t)
%% expoential map of D
I = speye(size(D));
tD = t*D;
tD2= t*t*D*D;

T = I + tD + 0.5*tD2;

end

