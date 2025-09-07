function P = assembleInnerPenalty1D(fem, grid, sigma)

    dPhi0 = fem.diffSpan([1,0], 1);
    dPhi1 = fem.diffSpan([0,1], 1);
    phi0  = fem.getSpan([1,0]);
    phi1  = fem.getSpan([0,1]);

    nDof = fem.locDof;
    NT   = length(grid) - 1;
    h    = grid(2:end) - grid(1:end-1);

    S1 = dPhi1' * phi1 + phi1' * dPhi1;
    M1 = phi1'  * phi1;

    S0 = - dPhi0' * phi0 - phi0' * dPhi0;
    M0 = phi0'  * phi0;

    C1 = - dPhi1' * phi0;
    C0 =   phi1' * dPhi0;
    J  = - phi1'  * phi0;

    hl    = h(1:NT-1);
    hr    = h(2:NT);
    hmean = 2 * hl .* hr ./ (hl + hr);
    a = 1 ./ (2 * hl);
    b = 1 ./ (2 * hr);
    w = sigma ./ (hmean.^2);

    ne = NT - 1;
    Id = speye(nDof);

    D_left_edge  = kron(spdiags(a, 0, ne, ne), S1) + kron(spdiags(w, 0, ne, ne), M1);
    D_right_edge = kron(spdiags(b, 0, ne, ne), S0) + kron(spdiags(w, 0, ne, ne), M0);
    OFF_edge     = kron(spdiags(a, 0, ne, ne), C1) + kron(spdiags(b, 0, ne, ne), C0) ...
                 + kron(spdiags(w, 0, ne, ne), J);

    A = sparse(1:ne, 1:ne, 1, NT, ne);
    B = sparse(2:NT, 1:ne, 1, NT, ne);

    LA = kron(A, Id);
    LB = kron(B, Id);

    P = sparse(NT*nDof, NT*nDof);
    P = P + LA * D_left_edge  * LA';
    P = P + LB * D_right_edge * LB';
    P = P + LA * OFF_edge * LB';
    P = P + (LA * OFF_edge * LB')';

end




%% 以下是手写的原始代码（无向量化优化），上面代码由 ChatGPT 根据手写代码优化而来
% function P = assembleInnerPenalty1D(fem, grid, sigma)
% 
%     dPhi0 = fem.diffSpan([1,0], 1);
%     dPhi1 = fem.diffSpan([0,1], 1);
%     phi0 = fem.getSpan([1,0]);
%     phi1 = fem.getSpan([0,1]);
% 
%     nDof = fem.locDof;
%     NT = length(grid) - 1;
%     h = grid(2:end) - grid(1:end-1);
%     P = sparse(NT*nDof, NT*nDof);
% 
%     for eid = 1 : NT-1
%         Ploc = zeros(nDof * 2);
%         hmean = (h(eid) + h(eid+1)) / 2;
%         %% u,v in K_left
%         dphi_mean = dPhi1 / (h(eid) * 2);
%         phi_jump = -phi1;
%         Ploc(1:nDof,1:nDof) = - dphi_mean' * phi_jump - phi_jump' * dphi_mean ...
%                               + (sigma/hmean) * (phi_jump' * phi_jump);
%         %% u,v in K_right
%         dphi_mean = dPhi0 / (h(eid+1) * 2);
%         phi_jump = phi0;
%         Ploc(nDof+1:end,nDof+1:end) = - dphi_mean' * phi_jump - phi_jump' * dphi_mean ...
%                                     + (sigma/hmean) * (phi_jump' * phi_jump);
%         %% u in K_left, v in K_right
%         du_mean = dPhi1 / (h(eid) * 2);
%         dv_mean = dPhi0 / (h(eid+1) * 2);
%         u_jump = -phi1;
%         v_jump = phi0;
%         Ploc(1:nDof,nDof+1:end) = - du_mean' * v_jump - u_jump' * dv_mean ...
%                                 + (sigma/hmean) * (u_jump' * v_jump);
%         %% u in K_right, v in K_left
%         Ploc(nDof+1:end,1:nDof) = Ploc(1:nDof,nDof+1:end)';
%         %% Assemble to global
%         idx = (eid-1)*nDof + (1 : 2*nDof);
%         P(idx, idx) = P(idx, idx) + Ploc;
%     end
% 
% end