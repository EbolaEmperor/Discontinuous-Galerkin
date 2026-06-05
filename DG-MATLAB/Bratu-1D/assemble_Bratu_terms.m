function [F, J] = assemble_Bratu_terms(fem, grid, u)
%ASSEMBLE_BRATU_TERMS Assemble nonlinear exp(u) terms for 1D Bratu equation.
%   F_i  = \int exp(u_h) * phi_i dx
%   J_ij = \int exp(u_h) * phi_i * phi_j dx

    NT = numel(grid) - 1;
    ndof = fem.nDof(NT);

    [quadL, w] = quadpts1_my(2 * fem.locDof);
    w = w(:)';
    nq = numel(w);
    phi = zeros(fem.locDof, nq);
    for i = 1:nq
        phi(:, i) = fem.getSpan(quadL(i, :))';
    end

    F = zeros(ndof, 1);
    J = sparse(ndof, ndof);
    for t = 1:NT
        vtx = grid(t:t+1);
        h = vtx(2) - vtx(1);
        idx = fem.dofMap(t);
        uhq = (u(idx).') * phi;
        eq = exp(uhq);
        F(idx) = F(idx) + sum(phi .* eq .* w, 2) * h;
        J(idx, idx) = J(idx, idx) + (phi .* eq .* w) * phi.' * h;
    end
end

