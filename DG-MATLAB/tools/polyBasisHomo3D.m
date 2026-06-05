function val = polyBasisHomo3D(k, p)
    number = numSplit3(k);
    val = (p(1).^number(1,:)) .* (p(2).^number(2,:)) .* (p(3).^number(3,:));
end