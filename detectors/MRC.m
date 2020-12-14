% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% Maximum Ratio Combining (MRC) detector
function [idxhat,bithat] = MRC(par,H,y)
  shat = H'*y;
  G = real(diag(H'*H));
  [~,idxhat] = min(abs(shat*ones(1,length(par.symbols))-G*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);
end
