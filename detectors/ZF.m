% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% unbiased Zero Forcing (ZF) detector
function [idxhat,bithat] = ZF(par,H,y)
  W = (H'*H)\(H');
  shat = W*y;
  G = real(diag(W*H));
  [~,idxhat] = min(abs(shat*ones(1,length(par.symbols))-G*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);
end
