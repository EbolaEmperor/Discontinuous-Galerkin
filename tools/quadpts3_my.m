function [lambda,weight] = quadpts3_my(order)

if order <= 5
    [lambda,weight] = quadpts3(order);
else 
    if order<=8
        [lambda,weight] = quadtetra(8);
    else
        [lambda,weight] = quadtetra(14);
    end
end

end