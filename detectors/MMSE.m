% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% unbiased Minimum Mean Squared-Error (MMSE) detector
function shat = MMSE(par,H,y,N0)
  W = (H'*H+(N0/par.Es)*eye(par.MT))\(H');
  shat = W*y;
  G = real(diag(W*H));  
  shat = shat./G;
end
