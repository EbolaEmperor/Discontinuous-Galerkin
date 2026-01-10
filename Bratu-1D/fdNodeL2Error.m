function e = fdNodeL2Error(u, x, h, u_exact)
%FDNODEL2ERROR Node-based L2 error on a uniform grid (via sqrt(h)*||.||_2).

    e = sqrt(h) * norm(u - u_exact(x));
end

