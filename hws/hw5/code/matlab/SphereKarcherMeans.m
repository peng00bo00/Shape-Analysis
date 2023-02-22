clc; clear;

N = 5;
Z = abs(randn(N,2)); Z = Z./vecnorm(Z,2,2);

M = obliquefactory(2, 1);
problem.M = M;
problem.cost = @(x,store)karchercost(x,store,Z);
problem.grad = @(x,store)karchergrad(x,store,Z);
problem.hess = @(x,v,store)karcherhess(x,v,store,Z);

x0 = [-1 -1]'/sqrt(2);
[xstar1, energy1, info1] = trustregions(problem, x0);
[xstar2, energy2, info2] = steepestdescent(problem, x0);

figure;
checkgradient(problem); disp(' ');
figure;
checkhessian(problem);

function [cost, store] = karchercost(x, store, Z)
    cost = sum(acos(Z*x).^2);
end

function [rgrad, store] = karchergrad(x, store, Z)
    %% Homework: Implement the Riemannian gradient as rgrad
    % euclidean gradient
    u = Z * x; d = acos(u);
    pdpx = -Z'./sqrt(1-u.*u)';
    egrad = 2*pdpx*d;
    
    %% projected to tangent plane
    rgrad = proj(x, egrad);
end

function [rhessv, store] = karcherhess(x, v, store, Z)
    %% Homework: Implement the Riemannian hessian as rhess
    rhess = 0;

    u = Z * x; d = acos(u);
    pdpx = -Z'./sqrt(1-u.*u)';
    egrad = 2*pdpx*d;

    ehess = Z' * diag(1./(1-u.*u)) * Z - Z' * diag(d.*u./((1-u.*u).^1.5)) * Z;
    ehess = 2 * ehess;

    rhess = proj(x, ehess) - x'*egrad;
    
    %% Multiply Riemannian Hessian with v
    rhessv = rhess*v;

end

function PX = proj(X, H)
%% project H to the tangent space of X
    PX = H - X * X' * H;
end