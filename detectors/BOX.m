% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% BOX detector
function [idxhat,bithat] = BOX(par,H,y)

  % -- initialization
  alpha = max(real(par.symbols));
  shat = H'*y;
  
  % -- apply a projected gradient descent
  for ii=1:par.BOX.iters
    shat = shat-par.BOX.tau*H'*(H*shat-y);
    shat = projinf(shat,alpha);
  end
  
  % -- compute outputs
  [~,idxhat] = min(abs(shat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);
  
end

% project into alpha infinity-norm ball
function sproj = projinf(s,alpha)
  sr = real(s);
  sr = max(min(sr,alpha),-alpha);  
  si = imag(s);
  si = max(min(si,alpha),-alpha);
  sproj = sr + 1i*si;
end
