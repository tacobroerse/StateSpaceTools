%% Create lower diagonal matrix using Householder transformations
% Version 1.0 October 2014 Thomas 'no comments' Frederikse

function[workmat] = householder(inputmat)
[pt,qt] = size(inputmat);
workmat = inputmat;
for i=1:pt
    reducemat = workmat(i:end,i:end);
    [~,q] = size(reducemat);
    x = reducemat(1,:);
    if (any(x(2:end))~=0)
        alpha = norm(x);
        e = zeros(1,q);
        e(1) = 1;
        u = x-alpha*e;
        v = u./norm(u);
        Q = eye(q) - 2*(v'*v);
        Qn = eye(qt);
        Qn(i:end,i:end) = Q;
        workmat = workmat*Qn;
    end
end

end