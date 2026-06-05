function [lambda,weight] = quadpts2_my(order)

if order <= 9
    [lambda,weight] = quadpts(order);
else 
    [lambda,weight] = quadtria(14);
end

end