classdef sinsin

properties
    u_exact
    grad_u_exact
    laplace_u_exact
end

methods
    function obj = sinsin(center)
        obj.u_exact = @(p) sin(pi*(p(:,1)-center)) .* sin(pi*(p(:,2)-center));

        obj.grad_u_exact = @(p) [cos(pi*(p(:,1)-center)) .* sin(pi*(p(:,2)-center)), ...
                                 sin(pi*(p(:,1)-center)) .* cos(pi*(p(:,2)-center))] * pi;

        obj.laplace_u_exact = @(p) (-2*pi*pi) * sin(pi*(p(:,1)-center)) .* sin(pi*(p(:,2)-center));
    end
end

end