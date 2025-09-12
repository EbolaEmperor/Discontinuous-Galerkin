classdef sin2sin2

properties
    u_exact
    grad_u_exact
    laplace_u_exact
    lap_lap_u_exact
end

methods
    function obj = sin2sin2(center)
        % 记 X = pi*(x-center), Y = pi*(y-center)
        
        % u_exact = sin^2 X * sin^2 Y
        obj.u_exact = @(p) ...
            (sin(pi*(p(:,1)-center)).^2) .* (sin(pi*(p(:,2)-center)).^2);

        % df/dx = 2*pi*sinX*cosX * sin^2Y
        % df/dy = 2*pi*sin^2X * sinY*cosY
        obj.grad_u_exact = @(p) ...
            [ 2*pi*sin(pi*(p(:,1)-center)).*cos(pi*(p(:,1)-center)) .* (sin(pi*(p(:,2)-center)).^2), ...
              2*pi*(sin(pi*(p(:,1)-center)).^2) .* sin(pi*(p(:,2)-center)).*cos(pi*(p(:,2)-center)) ];

        % Δf = π^2( cos(2X) + cos(2Y) - 2 cos(2X) cos(2Y) )
        obj.laplace_u_exact = @(p) ...
            (pi^2) * ( cos(2*pi*(p(:,1)-center)) + cos(2*pi*(p(:,2)-center)) ...
                     - 2*cos(2*pi*(p(:,1)-center)).*cos(2*pi*(p(:,2)-center)) );

        % Δ^2 f = -4π^4( cos(2X) + cos(2Y) - 4 cos(2X) cos(2Y) )
        obj.lap_lap_u_exact = @(p) ...
            (-4*pi^4) * ( cos(2*pi*(p(:,1)-center)) + cos(2*pi*(p(:,2)-center)) ...
                         - 4*cos(2*pi*(p(:,1)-center)).*cos(2*pi*(p(:,2)-center)) );
    end
end

end