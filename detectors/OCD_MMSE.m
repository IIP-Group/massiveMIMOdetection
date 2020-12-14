% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% Optimized Coordinate Descent (OCD) MMSE version
% -- Michael Wu, Chris Dick, Joseph R. Cavallaro, and Christoph Studer,
% -- "High-Throughput Data Detection for Massive MU-MIMO-OFDM Using
% -- Coordinate Descent," 
% -- IEEE Transactions on Circuits and Systems I: Regular Papers,
% -- vol. 63, no. 12, pp. 2357-2367, Dec. 2016.
function [idxhat,bithat] = OCD_MMSE(par,H,y,N0)

  % -- initialization
  alpha = N0/par.Es; % MMSE regularization; original code had a 0.5 factor

  % -- preprocessing
  dinv = zeros(par.MT,1);
  p = zeros(par.MT,1);
  for uu=1:par.MT
    normH2 = norm(H(:,uu),2)^2;
    dinv(uu,1) = 1/(normH2+alpha);
    p(uu,1) = dinv(uu)*normH2;
  end

  r = y;
  zold = zeros(par.MT,1);
  znew = zeros(par.MT,1);
  deltaz = zeros(par.MT,1);

  % -- OCD loop
  for iters=1:par.OCD_MMSE.iters
    for uu=1:par.MT
      znew(uu) = dinv(uu)*(H(:,uu)'*r)+p(uu)*zold(uu);
      deltaz(uu) = znew(uu)-zold(uu);
      r = r - H(:,uu)*deltaz(uu);
      zold(uu) = znew(uu);
    end
  end

  % -- compute outputs
  [~,idxhat] = min(abs(znew*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
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
