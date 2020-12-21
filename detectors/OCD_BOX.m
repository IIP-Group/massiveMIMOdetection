% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% Optimized Coordinate Descent (OCD) BOX version
% -- Michael Wu, Chris Dick, Joseph R. Cavallaro, and Christoph Studer,
% -- "High-Throughput Data Detection for Massive MU-MIMO-OFDM Using
% -- Coordinate Descent," 
% -- IEEE Transactions on Circuits and Systems I: Regular Papers,
% -- vol. 63, no. 12, pp. 2357-2367, Dec. 2016.
function shat = OCD_BOX(par,H,y)

  % -- initialization
  alpha = max(real(par.symbols));

  % -- preprocessing
  dinv = zeros(par.MT,1);
  p = zeros(par.MT,1);
  for uu=1:par.MT
    normH2 = norm(H(:,uu),2)^2;
    dinv(uu,1) = 1/normH2;
    p(uu,1) = dinv(uu)*normH2;
  end

  r = y;
  zold = zeros(par.MT,1);
  znew = zeros(par.MT,1);
  deltaz = zeros(par.MT,1);

  % -- OCD loop
  for iters=1:par.OCD_BOX.iters
    for uu=1:par.MT
      tmp = dinv(uu)*(H(:,uu)'*r)+p(uu)*zold(uu);
      znew(uu) = projinf(tmp,alpha);
      deltaz(uu) = znew(uu)-zold(uu);
      r = r - H(:,uu)*deltaz(uu);
      zold(uu) = znew(uu);
    end
  end

  shat = znew;

end

% project into alpha infinity-norm ball
function sproj = projinf(s,alpha)
  sr = real(s);
  sr = max(min(sr,alpha),-alpha);  
  si = imag(s);
  si = max(min(si,alpha),-alpha);
  sproj = sr + 1i*si;
end
