% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% Maximum Ratio Combining (MRC) detector
function shat = MRC(par,H,y)
  shat = H'*y;
  G = real(diag(H'*H));
  shat = shat./G;
end
