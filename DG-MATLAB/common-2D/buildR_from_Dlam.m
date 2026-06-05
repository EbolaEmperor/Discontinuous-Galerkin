function R = buildR_from_Dlam(Dlam)
    % Dlam: 2 x 3 x NT, 行分别是对 x,y 的导，列是 lam1,lam2,lam3
    NT = size(Dlam, 3);
    R  = zeros(3, 6, NT);

    a1 = squeeze(Dlam(1,1,:)).'; a2 = squeeze(Dlam(1,2,:)).'; a3 = squeeze(Dlam(1,3,:)).';
    b1 = squeeze(Dlam(2,1,:)).'; b2 = squeeze(Dlam(2,2,:)).'; b3 = squeeze(Dlam(2,3,:)).';

    % Hxx = a^T H a
    R(1,1,:) = reshape(a1.^2, 1,1,NT);
    R(1,2,:) = reshape(a2.^2, 1,1,NT);
    R(1,3,:) = reshape(a3.^2, 1,1,NT);
    R(1,4,:) = reshape(2*a1.*a2, 1,1,NT);
    R(1,5,:) = reshape(2*a1.*a3, 1,1,NT);
    R(1,6,:) = reshape(2*a2.*a3, 1,1,NT);

    % Hyy = b^T H b
    R(2,1,:) = reshape(b1.^2, 1,1,NT);
    R(2,2,:) = reshape(b2.^2, 1,1,NT);
    R(2,3,:) = reshape(b3.^2, 1,1,NT);
    R(2,4,:) = reshape(2*b1.*b2, 1,1,NT);
    R(2,5,:) = reshape(2*b1.*b3, 1,1,NT);
    R(2,6,:) = reshape(2*b2.*b3, 1,1,NT);

    % Hxy = a^T H b
    R(3,1,:) = reshape(a1.*b1, 1,1,NT);
    R(3,2,:) = reshape(a2.*b2, 1,1,NT);
    R(3,3,:) = reshape(a3.*b3, 1,1,NT);
    R(3,4,:) = reshape(a1.*b2 + a2.*b1, 1,1,NT);
    R(3,5,:) = reshape(a1.*b3 + a3.*b1, 1,1,NT);
    R(3,6,:) = reshape(a2.*b3 + a3.*b2, 1,1,NT);

    % Voigt 权重：把 Hxy 行乘上 sqrt(2)
    R(3,:,:) = R(3,:,:) * sqrt(2);
end