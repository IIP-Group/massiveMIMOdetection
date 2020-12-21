% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% detection via Triangular Approximate SEmidefinite Relaxation (TASER)
% -- Oscar Castañeda, Tom Goldstein, and Christoph Studer,
% -- "Data Detection in Large Multi-Antenna Wireless Systems via
% -- Approximate Semidefinite Relaxation," 
% -- IEEE Transactions on Circuits and Systems I: Regular Papers,
% -- vol. 63, no. 12, pp. 2334-2346, Dec. 2016.
function shat = TASER(par,H,y)

  switch par.mod
    case 'QPSK'
      % -- convert to real domain
      yR = [ real(y) ; imag(y) ];
      HR = [ real(H) -imag(H) ; imag(H) real(H) ];   
      % -- preprocessing for SDR  
      T = [HR'*HR , -HR'*yR ; -yR'*HR yR'*yR ];      
    case 'BPSK'
      % -- convert to real domain
      yR = [ real(y) ; imag(y) ];
      HR = [ real(H) ; imag(H) ];  
      % -- preprocessing for SDR  
      T = [HR'*HR , -HR'*yR ; -yR'*HR yR'*yR ];      
    otherwise
      error('modulation not supported')
  end
     
  DInv = diag(diag(T).^-.5);
  Ttilde = DInv*T*DInv;
  stepsize = par.TASER.alphaScale/norm(Ttilde,2);

  % -- use standard gradient on non-convex problem  
  gradf = @(L) 2*tril(L*Ttilde);
  proxg = @(L,t) prox_normalizer(L,diag(DInv).^-1);
  
  % Initialize Ltilde  
  Ltilde = diag(diag(DInv).^-1);
  
  % -- Fast Iterative Soft Thresholding [Beck & Tebouille, 2009]   
  for k = 1:par.TASER.iters
    Ltilde = proxg(Ltilde-stepsize*gradf(Ltilde)); % compute proxy    
  end  
  
  % -- post processing
  sRhat = sign(Ltilde(end,:))';  
  switch par.mod
    case 'QPSK'
      shat = sRhat(1:par.MT,1)+1i*sRhat(par.MT+1:end-1,1);
    case 'BPSK'  
      shat = sRhat(1:par.MT,1);
    otherwise
      error('modulation not supported')
  end
  
end

% normalize columns of Z to have norm equal to its corresponding scale
function Q = prox_normalizer(Z,scale)
  [N,~] = size(Z); 
  Q = Z.*(ones(N,1)*(1./sqrt(sum(abs(Z).^2,1)).*scale'));  
end
