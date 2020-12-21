% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% LAMA (Large MIMO Approximate Message Passing)
% Special thanks to Charles Jeon for this function
% -- Charles Jeon, Oscar Castañeda, and Christoph Studer
% -- "A 354 Mb/s 0.37 mm2 151 mW 32-User 256-QAM Near-MAP Soft-Input
% -- Soft-Output Massive MU-MIMO Data Detector in 28nm CMOS," 
% -- IEEE Solid-State Circuits Letters, vol. 2, no. 9, pp. 127-130, 
% -- Oct. 2019.
%
% For original algorithm, please see:
% -- Charles Jeon, Ramina Ghods, Arian Maleki, Christoph Studer
% -- "Optimality of large MIMO detection via approximate message passing,"
% -- IEEE International Symposium on Information Theory (ISIT),
% -- pp. 1227-1231, June 2015.
function shat = LAMA(par,H,y,N0)

  % transform MR x MT system to MT x MT by using Gram/matched filter
  % -- compute Gram matrix
  G = H'*H;
  % -- normalize Gram
  Gtilde = eye(par.MT) - diag(1./diag(G))*G;
  % -- compute column gains
  g = diag(G)/par.MR;
  % -- compute inverse col gain
  gtilde = 1./diag(G);
  % -- compute matched filter
  yMF = H'*y;
  % -- compute normalized matched filter
  yMFtilde = gtilde.*yMF;

  % Detector parameters
  % -- non-linear MMSE estimate
  shat   = zeros(par.MT,par.LAMA.iters+1);
  tau_s  = zeros(par.MT,par.LAMA.iters+1);
  tau_p  = zeros(1 ,    par.LAMA.iters+1);
  % -- signal mean/var estimate (modeled as s + Gaussian noise)
  z      = zeros(par.MT,par.LAMA.iters+1);
  tau_z = zeros(par.MT,par.LAMA.iters+1);
  % -- Onsager term
  v     = zeros(par.MT,par.LAMA.iters+1);

  % initialize estimates
  % -- assume input signal has variance par.Es
  tau_s(:,1) = par.Es * ones(par.MT,1);
  tau_p(1) = g' * tau_s(:,1);
  % -- first signal estimate is matched filter
  z(:,1) = yMFtilde;
  tau_z(:,1) = (tau_p(1) + N0)*gtilde;
  % -- damping parameters
  theta_tau_s = par.LAMA.theta_tau_s;
  theta_tau_z = par.LAMA.theta_tau_z;

  % loop
  for k=1:par.LAMA.iters    
    % -- Compute moments - mean/var
    [shat(:,k+1), tau_s(:,k+1)] = LAMA_MeanVar(par, z(:,k), tau_z(:,k));
    % -- damp second moment estimate
    tau_p(k+1) = theta_tau_s * ( g' * tau_s(:,k+1) ) ...
      + (1 - theta_tau_s) * tau_p(k);
    % -- Onsager term
    v(:,k+1) = tau_p(k+1) / (tau_p(k) + N0) * ( z(:,k) - shat(:,k) );
    % -- damp signal variance estimate
    tau_z(:,k+1) = theta_tau_z * (tau_p(k+1) + N0) * gtilde ...
      + (1-theta_tau_z) * tau_z(:,k);
    % -- signal update
    z(:,k+1) = yMFtilde + Gtilde * shat(:,k+1) + v(:,k+1);
  end

  % -- output signal
  shat = z(:,end);

end

% Special thanks to Charles Jeon for this function
function [F,G] = LAMA_MeanVar(par, z, sigma_sq)

  switch par.mod
    % for BPSK in complex-domain
    case 'BPSK'
        
      F = tanh(2*real(z)./(sigma_sq));
      G =  1 - abs(F).^2;        
      % compute estimate while considering numerical accuracy
      
    otherwise  
        
      % -- compute input matrices
      inputM = z*ones(1,length(par.symbols));
      
      % -- compute symbol matrices
      symbolM = ones(length(z),1)*par.symbols;
      sigma_sqM = sigma_sq*ones(1,length(par.symbols));
      
      % -- compute distance
      input_symM = abs(inputM - symbolM).^2;
      [~, index_vec] = min(input_symM,[],2);        
      symbolMinM = (par.symbols(index_vec).')*(ones(1,length(par.symbols)));
      expo = -(input_symM - abs(inputM - symbolMinM).^2);
      
      % -- compute pdf
      pdf = exp(expo./(sigma_sqM));
      pdf(isnan(pdf)) = 1;
      sumM = sum(pdf,2)*ones(1,length(par.symbols));
      
      % -- compute weight constant
      w = pdf./sumM;

      % -- compute mean (F)
      F = w*(par.symbols.');
      symbolM2 = ones(length(F),1)*par.symbols;
      FM = F*ones(1,length(par.symbols));
        
      % -- compute variance (G)
      G = sum(w.*(abs(symbolM2-FM).^2),2);
        
      if any(isnan(pdf))|any(isinf(pdf))
        warning('LAMA_MeanVar: NaN or inf OCCURRED');
      end

  end

end
