% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% SIMO lower bound
function shat = SIMO(par,H,y,s)
  z = y-H*s;
  shat = zeros(par.MT,1);
  for m=1:par.MT
    hm = H(:,m);
    yhat = z+hm*s(m,1);
    shat(m,1) = hm'*yhat/norm(hm,2)^2;    
  end  
end
