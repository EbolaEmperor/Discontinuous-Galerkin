function H = polyBasisHomoHess3D(k, p)
    number = numSplit3(k);
    a = number(1,:); b = number(2,:); c = number(3,:);

    p1 = p(1); p2 = p(2); p3 = p(3);

    H11 = a.*max(a-1,0) .* (p1.^max(a-2,0)) .* (p2.^b)          .* (p3.^c);
    H22 = b.*max(b-1,0) .* (p1.^a)          .* (p2.^max(b-2,0)) .* (p3.^c);
    H33 = c.*max(c-1,0) .* (p1.^a)          .* (p2.^b)          .* (p3.^max(c-2,0));

    H12 = a.*b .* (p1.^max(a-1,0)) .* (p2.^max(b-1,0)) .* (p3.^c);
    H13 = a.*c .* (p1.^max(a-1,0)) .* (p2.^b)          .* (p3.^max(c-1,0));
    H23 = b.*c .* (p1.^a)          .* (p2.^max(b-1,0)) .* (p3.^max(c-1,0));

    H = [H11; H22; H33; H12; H13; H23];
end