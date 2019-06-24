%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Accordance and Discordance %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to call the function, type con=compute(X,u,l) after defining X,u and l.
%X is the design matrix(time courses, regions), u,l the upper and lower 
%thresholds for Acc and Dis.
%eg. X = (std).*randn(10,20)+(mean); con=compute_con(X,1,1)

function [con] = compute_con(X,u,l)
N=size(X,1);% N: regions (in X, the rows represent the regions)
T=size(X,2);% T: time (in X, the columns represent the time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Normalization function %%%%%%%%%%%%%%%%%%%%%%%%%%
  function normalize(X)
        X = X - mean(X);
        X = X / std(X);
      for i=1:N
        X(i,:) = normalize(X(i,:)); %normalize each line of X matrix 
      end 
  end
%%%%%%%%%%%%%%%%Create zero matrices (Xu,Xl)%%%%%%%%%%%%%%%%%%
Xu=zeros(N,T);
Xl=zeros(N,T);
%%%%%%%%%%%%%%%%%%%%% Apply Thresholds %%%%%%%%%%%%%%%%%%%%%%%
Xu(X>=u)=1; %apply thresholds(u) on X and put ones(=1) in Xu for the corresponding positions 
Xl(X<=l)=-1;%apply theshold l on X and put -1 in Xl for the corresponding positions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% ACCORDANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Scav= 1 / T*(Xu *(Xu'));
Scdv= 1 / T*(Xl *(Xl'));
Sa=Scav + Scdv;
E=diag(Sa).^(-0.5); %energy
Acc=diag(E)* Sa * diag(E); %Calculate the Accordance matrix (Acc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% DISCORDANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sd = 1 / T*((Xu * (Xl')) + (Xl * (Xu')));
Dis = diag(E) * Sd * diag(E); %Calculate the Discordance matrix (Dis)

con=[];
con.accordance=Acc;
con.discordance=Dis;

end
 