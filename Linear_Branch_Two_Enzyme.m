function numerical_test_MM_robustness_branch1()
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
  %

 
  idx=find(1+Ps(:,2)-Ps(:,end).*Ps(:,2)>=0);
  Ps=Ps(idx,:);
  NN=size(Ps,1);
  clear idx
  %%
    Epsilon_rd=10.^(rand(NN,1).*5-10);
  Eo_max=25e-3/(2e-6); % maximum enzyme concentration in condensates
  PPs=[];
  ct=1;

  for i=1:NN   
      gamma=Ps(i,2);
      alpha=Ps(i,end);
      Eo=alpha*Ps(i,3);
      Ec=(1+gamma)*Ps(i,3);
      if Eo<=Eo_max/3 && Ec<=Eo_max/3 
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
   SR=RR./BR;
Res=[In' G theta RR' SR'];
save('FigS2A.mat','Res')
   figure(1)
   scatter(G,SR)
   
   %%
  Theta=[0.1 1 10];
  Res1=[];
  Res2=[];
  Res3=[];
  Res4=[];
  Res5=[];
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
      if 3*Eo<=Eo_max && 3*Ec<=Eo_max
         PPst=[PPst;Ps(i,:)];     


         [RR1(ct),BR1(ct),Ic1(ct),Io1(ct),In1(ct)]=ODE_eval2(Pss,1);
     
         ct=ct+1;
      end
   end
          Res1=[Res1 RR1' ];
       Res2=[Res2 RR1'./BR1'];
       Res3=[Res3 BR1'];
       Res4=[Res4 Ic1' Io1' In1'];
       clear RR1 BR1 Ic1 Io1 In1
  end
    gamma=PPst(:,2);
    alpha=PPst(:,8);
    %%
       G1=1+alpha.*gamma.*(10-1);
       Res=[Res1 Res2 gamma.*alpha];
       save('Fig3A.mat','Res')
   
   
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






%% three enzyme linear+branch pathway
x0=zeros(4,1);
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
Yn=yn(end,3);
In=yn(end,1);
Bn=yn(end,end);
x0=zeros(8,1);
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
Ass=y(end,5)/(1+gamma)+(gamma/(1+gamma))*(y(end,6));
YSS=Ass/Yn;

Ic=y(end,1);
Io=y(end,2);
Bss=y(end,end-1)/(1+gamma)+(gamma/(1+gamma))*(y(end,end));
BSS=Bss/Bn;



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


 


%% three enzyme linear+branch system: limited upstream metabolite
Sc=x(1);
So=x(2);
Ic=x(3);
Io=x(4);
Ac=x(5);
Ao=x(6);
Bc=x(7);
Bo=x(8);
E1o=E*alpha;
E1c=E*(1+gamma - alpha*gamma);
E2o=E*alpha;
E2c=E*(1+gamma - alpha*gamma);
E3c=E*(1+gamma);
dScdt=k0             - k*E1c*Sc/(1+Sc)                  -      beta*Sc - eta*(Sc-So)*gamma;
dSodt=               - theta*k*E1o*So/(1+So)            -      beta*So - eta*(So-Sc);

dIcdt=k*E1c*Sc/(1+Sc)       - k*E2c*Ic/(1+Ic)       - k*E3c*Ic/(1+Ic) - beta*Ic - eta*(Ic-Io)*gamma;
dIodt=theta*k*E1o*So/(1+So) - theta*k*E2o*Io/(1+Io)            - beta*Io - eta*(Io-Ic);

dAcdt=                 k*E2c*Ic/(1+Ic)                  -      Ac - eta*(Ac-Ao)*gamma;
dAodt=                 theta*k*E2o*Io/(1+Io)            -      Ao - eta*(Ao-Ac);

dBcdt=                 k*E3c*Ic/(1+Ic)                  -      Bc - eta*(Bc-Bo)*gamma;
dBodt=                                           -      Bo - eta*(Bo-Bc);
J=[dScdt;dSodt;dIcdt;dIodt;dAcdt;dAodt;dBcdt;dBodt];
%
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



%% three enzyme linear+branch system: 
Sn=x(1);
In=x(2);
An=x(3);
Bn=x(4);
dSndt=k0-k*E*Sn/(1+Sn)-beta*Sn;
dIndt=   k*E*Sn/(1+Sn)-2*k*E*In/(1+In)-beta*In;
dAndt=   k*E*In/(1+In)-An;
dBndt=   k*E*In/(1+In)-Bn;
J=[dSndt;dIndt;dAndt;dBndt];

end
