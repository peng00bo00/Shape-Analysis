addpath('utils/');

meshfile = 'meshes/166.off';
% meshfile = 'meshes/moomoo.off';
[X,T] = readOff(meshfile);
nv = size(X,1);
nt = size(T,1);

A = surfaceArea(X,T);
fprintf('The surface area of %s is %f\n', meshfile, A);

% Sanity checks: Laplacian is symmetric and positive definite
L = cotLaplacian(X,T);
[~,p] = chol(L);
fprintf('\nIf %d is 0, then L is PSD\n', p);
fprintf('Symmetry: %d\n', norm(L-L','fro'));

% Divided differences approximation
%% it may take a long time to compute numerical gradient
% G = gradApprox(X,T);
% fprintf('Difference between gradient and cotan Laplacian %f\n', norm(-.5*L*X-G,'fro'));

% Barycentric area
fprintf('Difference between sum of barycentric area and surface area %f\n', sum(barycentricArea(X,T))-A);

H = meanCurvature(X,T);
showDescriptor(X,T,H);

%% Mean curvature flow
maxiters=1000;
Xt = X;
% ADD CODE FOR EXPLICIT INTEGRATOR HERE %%%%
tau = 0.000005;
saveas(gcf,['explicit_integrator_step-0.png']);
for t=1:maxiters
    t
    cotL =-cotLaplacian(Xt,T);
    Xt = Xt - tau*cotL*Xt./barycentricArea(Xt,T);

    if mod(t, 10) == 0
        H = meanCurvature(Xt,T);
        showDescriptor(Xt,T,H);
        saveas(gcf,['explicit_integrator_step-' int2str(t) '.png']);
    end
end
% END HOMEWORK ASSIGNMENT %%%%
% Uncomment to display mesh at the end
H = meanCurvature(Xt,T);
showDescriptor(Xt,T,H);

Xt = X;
% ADD CODE FOR IMPLICIT INTEGRATOR HERE %%%%
maxiters=1000;
tau = 0.00001;
for t=1:maxiters
    t
    cotL =-cotLaplacian(Xt,T);
    A = speye(nv)+tau*cotL./barycentricArea(Xt,T);
    Xt = A \ Xt;

    if mod(t, 10) == 0
        H = meanCurvature(Xt,T);
        showDescriptor(Xt,T,H);
        saveas(gcf,['implicit_integrator_step-' int2str(t) '.png']);
    end
end
% END HOMEWORK ASSIGNMENT %%%%
% Uncomment to display mesh at the end
H = meanCurvature(Xt,T);
showDescriptor(Xt,T,H);

%% Non-singular Mean Curvature Flow
Xt = X;
maxiters=1000;
tau = 0.00001;
cotL =-cotLaplacian(Xt,T);

for t=1:maxiters
    t
    A = speye(nv)+tau*cotL./barycentricArea(Xt,T);
    Xt = A \ Xt;

    if mod(t, 10) == 0
        H = meanCurvature(Xt,T);
        showDescriptor(Xt,T,H);
        saveas(gcf,['non-sigular_implicit_integrator_step-' int2str(t) '.png']);
    end
end
H = meanCurvature(Xt,T);
showDescriptor(Xt,T,H);

%% Function definitions
% ADD CODE TO COMPUTE SURFACE AREA HERE %%%%%%%%%%
function [A] = surfaceArea(X,T)
    A = 0;

    for i=1:size(T, 1)
        A = A + triangleArea(X(T(i,:),:));
    end

end
% END HOMEWORK ASSIGNMENT %%%%%%%%%%%

% ADD CODE TO COMPUTE COTANGENT LAPLACIAN HERE %%%%%%%%%
function [L] = cotLaplacian(X,T)
    nv = size(X,1);
    nt = size(T,1);

    %% indices and values
    I = zeros(nt*12,1);
    J = zeros(nt*12,1);
    V = zeros(nt*12,1);

    for t = 1:nt
        %% vertices
        X1 = X(T(t,1),:);
        X2 = X(T(t,2),:);
        X3 = X(T(t,3),:);

        %% edges
        e12= X2 - X1;
        e23= X3 - X2;
        e31= X1 - X3;

        %% cot
        cot1 = cotTheta(e12,-e31);
        cot2 = cotTheta(e23,-e12);
        cot3 = cotTheta(e31,-e23);

        %% check http://rodolphe-vaillant.fr/entry/20/compute-harmonic-weights-on-a-triangular-mesh for more details
        %% also note that we'll accumulate values on each vertex
        
        %% e23
        I(t*12-11) = T(t,2); J(t*12-11) = T(t,3); V(t*12-11) = cot1;
        I(t*12-10) = T(t,3); J(t*12-10) = T(t,2); V(t*12-10) = cot1;
        I(t*12-9)  = T(t,2); J(t*12-9)  = T(t,2); V(t*12-9)  =-cot1;
        I(t*12-8)  = T(t,3); J(t*12-8)  = T(t,3); V(t*12-8)  =-cot1;

%         L(T(t,2), T(t,3)) = L(T(t,2), T(t,3)) + cot1;
%         L(T(t,3), T(t,2)) = L(T(t,3), T(t,2)) + cot1;
%         L(T(t,2), T(t,2)) = L(T(t,2), T(t,2)) - cot1;
%         L(T(t,3), T(t,3)) = L(T(t,3), T(t,3)) - cot1;
        
        %% e31
        I(t*12-7)  = T(t,3); J(t*12-7)  = T(t,1); V(t*12-7)  = cot2;
        I(t*12-6)  = T(t,1); J(t*12-6)  = T(t,3); V(t*12-6)  = cot2;
        I(t*12-5)  = T(t,3); J(t*12-5)  = T(t,3); V(t*12-5)  =-cot2;
        I(t*12-4)  = T(t,1); J(t*12-4)  = T(t,1); V(t*12-4)  =-cot2;

%         L(T(t,3), T(t,1)) = L(T(t,3), T(t,1)) + cot2;
%         L(T(t,1), T(t,3)) = L(T(t,1), T(t,3)) + cot2;
%         L(T(t,3), T(t,3)) = L(T(t,3), T(t,3)) - cot2;
%         L(T(t,1), T(t,1)) = L(T(t,1), T(t,1)) - cot2;

        %% e12
        I(t*12-3)  = T(t,1); J(t*12-3)  = T(t,2); V(t*12-3)  = cot3;
        I(t*12-2)  = T(t,2); J(t*12-2)  = T(t,1); V(t*12-2)  = cot3;
        I(t*12-1)  = T(t,1); J(t*12-1)  = T(t,1); V(t*12-1)  =-cot3;
        I(t*12)    = T(t,2); J(t*12)    = T(t,2); V(t*12)    =-cot3;

%         L(T(t,1), T(t,2)) = L(T(t,1), T(t,2)) + cot3;
%         L(T(t,2), T(t,1)) = L(T(t,2), T(t,1)) + cot3;
%         L(T(t,1), T(t,1)) = L(T(t,1), T(t,1)) - cot3;
%         L(T(t,2), T(t,2)) = L(T(t,2), T(t,2)) - cot3;

    end

    %% build cotLaplacian
    L = sparse(I,J,V,nv,nv);
end
% END HOMEWORK ASSIGNMENT %%%%%%%%%%%

% ADD CODE TO COMPUTE DIVIDED DIFFERENCES APPROXIMATION HERE %%%%%
function [G] = gradApprox(X,T)
    nv = size(X,1);
    G = zeros(nv,3);

    %% step size
    h = 0.001;

    for v = 1:nv
        %% x,y,z axis
        for i = 1:3
            Xplus = X;
            Xplus(v, i)= Xplus(v, i) + h;

            Xminus= X;
            Xminus(v,i)= Xminus(v,i) - h;

            G(v, i) = (surfaceArea(Xplus, T) - surfaceArea(Xminus, T)) / (2*h);
        end
    end

%     nv = size(X,1);
%     nt = size(T,1);
%     G = zeros(nv,3);
%     h = 0.001;
%     for axis = 1:3
%         % the change of area w.r.t each point is local
%         for t = 1:nt
%             v1 = X(T(t,1),:); v2 = X(T(t,2),:); v3 = X(T(t,3),:);
%             e = zeros(1,3);
%             e(axis) = 1;
%             G(T(t,1),axis) = G(T(t,1),axis)+(tarea(v1+e*h,v2,v3)-tarea(v1-e*h,v2,v3))/(2*h);
%             G(T(t,2),axis) = G(T(t,2),axis)+(tarea(v1,v2+e*h,v3)-tarea(v1,v2-e*h,v3))/(2*h);
%             G(T(t,3),axis) = G(T(t,3),axis)+(tarea(v1,v2,v3+e*h)-tarea(v1,v2,v3-e*h))/(2*h);
%         end
%     end
%     function tarea = tarea(X1,X2,X3)
%         tarea = norm(cross(X1-X2,X3-X2))/2;
%     end
end
% END HOMEWORK ASSIGNMENT %%%%%%%%%%%%%%%

% ADD CODE TO COMPUTE THE BARYCENTRIC AREA VECTOR HERE %%%%%%%%%
function [M] = barycentricArea(X,T)
    nv = size(X,1);
    nt = size(T,1);
    M = zeros(nv,1);

    for t = 1:nt    
        %% 1/3 of triangle area
        A3= triangleArea(X(T(t,:),:))/3;

        M(T(t,1)) = M(T(t,1)) + A3;
        M(T(t,2)) = M(T(t,2)) + A3;
        M(T(t,3)) = M(T(t,3)) + A3;
    end

end
% END HOMEWORK ASSIGNMENT %%%%%%%%%%%%%%

% ADD CODE TO COMPUTE POINTWISE MEAN CURVATURE HERE %%%%%%%%%%%
function [H] = meanCurvature(X,T)
    nv = size(X,1);
    H = zeros(nv,1);

    Hn= 0.5 * cotLaplacian(X,T) * X;
    H = sqrt(sum(Hn.*Hn, 2)) ./ barycentricArea(X,T);
end
% END HOMEWORK ASSIGNMENT %%%%%%%%%%%%

%% A HELPER FUNCTION TO CALCULATE TRIANGLE AREA
function [A] = triangleArea(T)
    X1 = T(1,:);
    X2 = T(2,:);
    X3 = T(3,:);

    A = 0.5*norm(cross(X2-X1, X3-X1));
end

%% A HELPER FUNCTION TO CALCULATE COT(THETA)
function [c] = cotTheta(v, w)
    c = dot(v, w) / norm(cross(v, w));
end