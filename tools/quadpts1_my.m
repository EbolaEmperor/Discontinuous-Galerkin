function [x, w] = quadpts1_my(n)
    i = (1:n-1)';
    b = i ./ sqrt(4*i.^2 - 1);
    J = diag(b,1) + diag(b,-1);
    [V,D] = eig(J);
    [x,perm] = sort(diag(D));
    V = V(:,perm);
    w = 2*(V(1,:).^2)';
    x = [(x+1)/2, (1-x)/2];
    w = w/2;
end