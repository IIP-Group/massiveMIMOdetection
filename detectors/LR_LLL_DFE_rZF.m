% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% detection via LLL-Lattice-Reduction with decision feedback and ZF remap
% -- Christoph Windpassinger and Robert F. H. Fischer,
% -- "Low-Complexity Near-Maximum-Likelihood Detection and Precoding for
% -- MIMO Systems using Lattice Reduction," Proceedings of the IEEE 
% -- Information Theory Workshop, Mar. 2003, pp. 345-348
function [idxhat,bithat] = LR_LLL_DFE_rZF(par,H,y)

  % -- initialization
  sEst = zeros(par.MT,1);
    
  % -- convert to lattice: s = alpha*(d+v)
  switch par.mod
    case 'BPSK'
      v = 0.5;
    case {'QPSK','16QAM','64QAM'}
      v = 0.5+1i*0.5;
    otherwise
      error('modulation not supported')
  end
  alpha = 2;
  d = par.symbols/alpha-v; % candidates in d-domain
  y_tilde = y-alpha*H*ones(par.MT,1)*v;
  H_tilde = H*alpha;
   
  % -- complex valued shifted and scaled constellations
  [Q,R,T] = LLL(par,H_tilde);
  y_hat = Q'*y_tilde;
  
  % -- DFE stages
  for Level=par.MT:-1:1
       
    % -- decision feedback (DF)
    DF = R(Level,Level+1:end)*sEst(Level+1:end);
    y_slice = y_hat(Level) - DF;
    
    %-- relaxed detection
    sEst(Level,1) = round(y_slice/R(Level,Level));
       
  end

  % -- ZF remap to original symbol
  sDFE = T*sEst;
  [~,idxhat] = min(abs(sDFE*ones(1,length(par.symbols))-ones(par.MT,1)*d).^2,[],2);
  bithat = par.bits(idxhat,:);  
  
end

function [Q,R,T,P]= LLL(par,H)
    
  [Q,R] = qr(H);
  P = eye(par.MT);
   
  T = eye(par.MT);  
  
  m = 2;
  while m <= par.MT
    % -- size reduction
    for n=m-1:-1:1
      mu = round(R(n,m)/R(n,n));
      R(1:n,m) = R(1:n,m)-mu*R(1:n,n);
      T(:,m) = T(:,m)-mu*T(:,n); 
    end    
    % -- L3 LR criterion
    if 0.75*abs(R(m-1,m-1))^2 > abs(R(m,m))^2 + abs(R(m-1,m))^2      
      % -- swap
      tmp = R(:,m-1); R(:,m-1) = R(:,m); R(:,m) = tmp;            
      tmp = T(:,m-1); T(:,m-1) = T(:,m); T(:,m) = tmp;
      
      % -- Givens rotation 
      Norm = sqrt(abs(R(m-1:m,m-1)'*R(m-1:m,m-1)));
      c = R(m-1,m-1)/Norm; s = R(m,m-1)/Norm;
      G = [c' s'; -s c];
      R(m-1:m,m-1:par.MT) = G*R(m-1:m,m-1:par.MT);
      Q(:,m-1:m) = Q(:,m-1:m)*G';
      
      m = max(m-1,2);      
    else
      m = m+1;
    end
  end

end
