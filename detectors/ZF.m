% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% unbiased Zero Forcing (ZF) detector
function shat = ZF(H,y)
  W = (H'*H)\(H');
  shat = W*y;
  G = real(diag(W*H));
  shat = shat./G;
end
