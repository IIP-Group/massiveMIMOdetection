% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% K-Best detector
function [idxhat,bithat] = KBEST(par,H,y)

  % -- preprocessing
  [Q,R] = qr(H);
  y_hat = Q'*y;

  % -- Initialize Partial Euclidean Distance (PED) with last TX symbol
  PED_list=abs(par.symbols*R(par.MT,par.MT) - y_hat(par.MT)).^2;
  [PED_list,idx]=sort(PED_list);
  s=par.symbols(:,idx);
  % -- take the K-best
  s=s(:,1:min(par.KBEST.K,length(PED_list)));
  Kbest_PED_list=PED_list(1:min(par.KBEST.K,length(PED_list)));

  % -- for each TX symbol
  for Level=par.MT-1:-1:1
    PED_list=[];
    % -- obtain the cumulative Euclidean distance considering the K-best
    %    previous nodes
    for k=1:length(Kbest_PED_list)
      tmp=Kbest_PED_list(k)+abs(par.symbols*R(Level,Level)-y_hat(Level) + ...
          R(Level,Level+1:par.MT)*s(:,k)).^2;
      PED_list=[PED_list,tmp];
    end
    % -- sort in ascending order
    s=[kron(ones(1,length(Kbest_PED_list)),par.symbols); ...
       kron(s,ones(1,length(par.symbols)))];
    [PED_list,idx]=sort(PED_list);
    s=s(:,idx);
    % take the K-best
    s=s(:,1:min(par.KBEST.K,length(PED_list)));
    Kbest_PED_list=PED_list(1:min(par.KBEST.K,length(PED_list)));
  end
  % -- take the best
  s=s(:,1);
  
  % -- compute outputs
  idxhat=zeros(par.MT,1);
  for i=[1:par.MT]
    idxhat(i,1)= find(s(i)==par.symbols);
  end  
  bithat = par.bits(idxhat,:);
  
end
