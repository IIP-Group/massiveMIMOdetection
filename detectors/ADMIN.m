% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% ADMM-based Infinity-Norm detection (ADMIN)
% -- Shariar Shahabuddin, Markku Juntti, and Christoph Studer,
% -- "ADMM-based Infinity Norm Detection for Large MU-MIMO: Algorithm and
% -- VLSI Architecture," 
% -- IEEE International Symposium on Circuits, Systems (ISCAS),
% -- May 2017.
function [idxhat,bithat] = ADMIN(par,H,y,N0)

  % -- initialization
  beta = N0/par.Es*par.ADMIN.betaScale; 
  G = H'*H + beta*eye(par.MT);
  [L,D] = ldl(G);
  Dinv = diag(1./diag(D));
  yMF = H'*y;
  zhat = zeros(par.MT,1);
  lambda = zeros(par.MT,1);
  alpha = max(real(par.symbols));
  
  % -- main loop
  for k = 1:par.ADMIN.iters
      shat = L'\(Dinv*(L\(yMF+beta*(zhat-lambda))));
      zhat = projinf(shat+lambda,alpha);
      lambda = lambda-par.ADMIN.gamma*(zhat-shat);
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
