% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% detection via TASER followed by randomization (TASER_R)
% -- Oscar Castañeda, Tom Goldstein, and Christoph Studer,
% -- "Data Detection in Large Multi-Antenna Wireless Systems via
% -- Approximate Semidefinite Relaxation," 
% -- IEEE Transactions on Circuits and Systems I: Regular Papers,
% -- vol. 63, no. 12, pp. 2334-2346, Dec. 2016.
function shat = TASER_R(par,H,y)

  switch par.mod
    case 'QPSK'
      % -- convert to real domain
      yR = [ real(y) ; imag(y) ];
      HR = [ real(H) -imag(H) ; imag(H) real(H) ];   
      % -- preprocessing for SDR  
      T = [HR'*HR , -HR'*yR ; -yR'*HR yR'*yR ];
      N = 2*par.MT+1;
    case 'BPSK'
      % -- convert to real domain
      yR = [ real(y) ; imag(y) ];
      HR = [ real(H) ; imag(H) ];  
      % -- preprocessing for SDR  
      T = [HR'*HR , -HR'*yR ; -yR'*HR yR'*yR ];
      N = par.MT+1;
    otherwise
      error('modulation not supported')
  end
     
  DInv = diag(diag(T).^-.5);
  Ttilde = DInv*T*DInv;
  stepsize = par.TASER_R.alphaScale/norm(Ttilde,2);

  % -- use standard gradient on non-convex problem  
  gradf = @(L) 2*tril(L*Ttilde);
  proxg = @(L,t) prox_normalizer(L,diag(DInv).^-1);
  
  % Initialize Ltilde  
  Ltilde = diag(diag(DInv).^-1);
  
  % -- Fast Iterative Soft Thresholding [Beck & Tebouille, 2009]   
  for k = 1:par.TASER_R.iters
    Ltilde = proxg(Ltilde-stepsize*gradf(Ltilde)); % compute proxy    
  end  
  
  % -- post processing with randomization
  z = randn(N,par.TASER_R.L);
  %z = sign(randn(N,par.TASER_R.L));
  %z = 2*randi(2,N,par.TASER_R.L)-3;
  LtildeT = Ltilde.';
  z = LtildeT*z;
  % we need last entry to be 1 after taking the sign
  z = z.*(ones(N,1)*sign(z(end,:)));
  z = sign(z);
  z_cost = diag(z'*T*z);
  [~,z_idx] = min(z_cost);
  
  sRhat = z(:,z_idx);  
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
