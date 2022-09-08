function y = mvu(x, k)

n = size(x, 2);
y = x;

%%% PROBLEM 4 - YOUR CODE HERE %%%

%% Gram matrix X
x = x - mean(x, 2);
X = x' * x;

%% distance matrix D
DX = diag(X) + diag(X)' - 2*X;

%% nearest neighbors
[~, N] = sort(DX, 2);
N = N(:,1:k+1);

%% neighbor graph
E = zeros(n*k*(k+1)/2, 2);

row = 1;
for v=1:n
    for ii=1:k+1
        i = N(v, ii);
        for jj=ii+1:k+1
            j=N(v, jj);
            
            E(row,1) = min(i, j);
            E(row,2) = max(i, j);

            row = row+1;
        end
    end
end

E = unique(E,'rows');

%% use CVX for solving
cvx_precision high
cvx_begin
    variable Y(n,n) semidefinite;

    maximize(trace(Y));

    subject to
        Y * ones(n, 1) == 0;

        for e=1:size(E, 1)
            i = E(e, 1);
            j = E(e, 2);

            Y(i,i)+Y(j,j)-2*Y(i,j) == DX(i,j);
        end

% %% KNN
%         for i=1:n
%             for jj=2:k+1
%                 j = N(i, jj);
%                 Y(i,i)+Y(j,j)-2*Y(i,j) == DX(i,j);
%             end
%         end
% cvx_end

[V,D] = eig(Y);
D = diag(D);
[D, ind] = sort(D, 'descend');
V = V(:, ind);

y = diag(sqrt(D)) * V';

%%% END HOMEWORK PROBLEM
    
end