function Fig1()
%P=  [ eta    gamma      E     theta   k    beta    k0     alpha   ]
P0=  [ 1e8     1e-2      0.02    1     1e5    1    2000      5];


Pmax=[1e9    1e-1       2       10    1e7    10     2e5     100];
Pmin=[1e5    1e-4     0.0002    0.1   1e3    0.1    20       1.3]; 


  NN=1000;
  p = sobolset(length(P0),'Skip',1e3,'Leap',1e2);
  p = scramble(p,'MatousekAffineOwen');
  Ps = p(1:3:3*NN,:);
  Ps=10.^(Ps.*repmat(log10(Pmax)-log10(Pmin),NN,1)+repmat(log10(Pmin),NN,1));
  clear p
 
  idx=find(1+Ps(:,2)-Ps(:,end).*Ps(:,2)>=0);
  Ps=Ps(idx,:);
  NN=size(Ps,1);
  clear idx
  %%
    
  Eo_max=25e-3/(2e-6); % maximum enzyme concentration in condensates
  PPs=[];
  ct=1;

  for i=1:NN   
      gamma=Ps(i,2);
      alpha=Ps(i,end);
      Eo=alpha*Ps(i,3);
      Ec=(1+gamma)*Ps(i,3);
      if Eo<=Eo_max && Ec<=Eo_max
         PPs=[PPs;Ps(i,:)];     
         [RR(ct),BR(ct),Ic(ct),Io(ct),In(ct)]=ODE_eval2(Ps(i,:),1);
        ct=ct+1;
      end
  end
     figure(1)

  eta=PPs(:,1);
   gamma=PPs(:,2);
   E=PPs(:,3);
   theta=PPs(:,4);
   k=PPs(:,5);
   beta=PPs(:,6);
   k0=PPs(:,7);
   alpha=PPs(:,8);
   Ek=E.*k;


   G=1+alpha.*gamma.*(theta-1);
   

   figure(1)
   scatter(G,RR)

 Res=[In' G theta RR'];
 save('Fig1b.mat','Res')
 %%
 Theta=[0.1 1 10];
 Res1=[];
Res4=[];
 for j=1:length(Theta)
    ct=1;
   PPst=[];
   for i=1:NN  
       Pss=P0;
    
      gamma=Ps(i,2);
      alpha=Ps(i,end);
      Eo=alpha*Ps(i,3);
      Ec=(1+gamma)*Ps(i,3);
         Pss(2)=gamma;
         Pss(end)=alpha;
         Pss(4)=Theta(j);
         Pss(7)=2e4;
      if Eo<=Eo_max && Ec<=Eo_max
         PPst=[PPst;Ps(i,:)];     

         [RR1(ct),BR1(ct),Ic1(ct),Io1(ct),In1(ct)]=ODE_eval2(Pss,1);
    
         ct=ct+1;
      end
   end
    gamma=PPst(:,2);
    alpha=PPst(:,8);
       G1=1+alpha.*gamma.*(Theta(j)-1);
       Res1=[Res1 RR1'];
      
       clear RR1 Ic1 Io1 In1
 end    
 clear Res
       Res=[Res1 alpha.*gamma];
       save('Fig1d.mat','Res')
   
end
function [YSS,BSS,Ic,Io,In]=ODE_eval2(P,Epsilon_rd)

eta=P(1);
gamma=P(2);
E=P(3);
theta=P(4);
k=P(5);
beta=P(6);
k0=P(7);
alpha=P(8);

 %% one enzyme linear pathway
x0=[0;0];
opts = odeset('RelTol',1e-3,'AbsTol',1e-3,'NonNegative',[1:length(x0)]);
flag=false;
while ~flag
    [~,yn]=ode15s(@(t,x) ODE_eval1_uncompart(t,x,P,Epsilon_rd),[0 100],x0,opts);
    dyn=ODE_eval1_uncompart(0,yn(end,:),P,Epsilon_rd);
    dyn_nor=abs(dyn./yn(end-5:end,:)');
    idx=find(dyn_nor<1e-2);
    if ~isempty(idx)
      flag=true;
    else
        x0=yn(end,:)';
    end
end
Yn=yn(end,2);
In=yn(end,1);
x0=[0;0;0;0];
opts = odeset('RelTol',1e-3,'AbsTol',1e-3,'NonNegative',[1:length(x0)]);
flag=false;
while ~flag
    [~,y]=ode15s(@(t,x) ODE_eval1(t,x,P,Epsilon_rd),[0 100],x0,opts);
    dy=ODE_eval1(0,y(end,:),P,Epsilon_rd);
    dy_nor=abs(dy./y(end-5:end,:)');
    idx=find(dy_nor<1e-2);
    if ~isempty(idx)
       flag=true;
    else
        x0=y(end,:)';
    end
end
Ass=y(end,3)/(1+gamma)+(gamma/(1+gamma))*(y(end,4));
YSS=Ass/Yn;

Ic=y(end,1);
Io=y(end,2);
BSS=0;


 
end
function J=ODE_eval1(t,x,P,Epsilon_rd)

%P=  [ eta    gamma      E     theta   k    beta    k0     alpha ]
eta=P(1);
gamma=P(2);
E=P(3);
theta=P(4);
k=P(5);
beta=P(6);
k0=P(7);
alpha=P(8);




%% One enzyme linear pathway
Ic=x(1);
Io=x(2);
Ac=x(3);
Ao=x(4);
Eo=E*alpha;
Ec=E*(1+gamma - alpha*gamma);
% % % dIcdt=k0 -      k*Ec*Ic/(1+Ic) - beta*Ic - eta*(Epsilon_rd*Ic-Io);
% % % dIodt=   -theta*k*Eo*Io/(1+Io) - beta*Io - eta*(Io-Epsilon_rd*Ic)/gamma;
% % % 
% % % dAcdt=          k*Ec*Ic/(1+Ic) -      Ac - eta*(Ac-Ao);
% % % dAodt=    theta*k*Eo*Io/(1+Io) -      Ao - eta*(Ao-Ac)/gamma;


dIcdt=k0 -      k*Ec*Ic/(1+Ic) - beta*Ic - eta*(Epsilon_rd*Ic-Io)*gamma;
dIodt=   -theta*k*Eo*Io/(1+Io) - beta*Io - eta*(Io-Epsilon_rd*Ic);

dAcdt=          k*Ec*Ic/(1+Ic) -      Ac - eta*(Ac-Ao)*gamma;
dAodt=    theta*k*Eo*Io/(1+Io) -      Ao - eta*(Ao-Ac);

J=[dIcdt;dIodt;dAcdt;dAodt];


end

function J=ODE_eval1_uncompart(t,x,P,Epsilon_rd)

%P=  [ eta    gamma      E     theta   k    beta    k0     alpha ]
eta=P(1);
gamma=P(2);
E=P(3);
theta=P(4);
k=P(5);
beta=P(6);
k0=P(7);
alpha=P(8);
epsilon=Epsilon_rd;

%% One enzyme linear pathway
In=x(1);
An=x(2);


dIndt=k0-k*E*In/(1+In)-beta*In;
dAndt=   k*E*In/(1+In)-An;
J=[dIndt;dAndt];


end
