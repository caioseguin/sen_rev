function A = sen_rev_asy_ttest(X)

% Function to compute communication asymmetry in X, a NxNxK communication efficiency matrix
% A is the resulting t-statistic matrix of pairwise communication asymmetries
% df is the degree of freedom of the t-tests performed
%
% Caio Seguin Feb 2019

    [n, m, k] = size(X);
    
    assert(k > 1, 'X must be a 3D matrix');
    assert(n == m, 'X must be a square matrix');
    
    A = zeros(n);
    
    for i = 1:n
        for j = 1:n
            if i ~= j
                x = squeeze(X(i,j,:) - X(j,i,:));
                A(i,j) = mean(x)/(std(x)/(sqrt(k)));
            end
        end
    end
    
    A(isnan(A)) = 0;

end