
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculates the equilibrium probability that a bp is open for
%% a homopolymer with a defect bp. The Poland algorithm is used,
%% and compared to analytics. The first and the last of the bps 
%% are assumed to be clamped.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
format long;


% ------- Input parameters -------

% Hydrogen bond energies
alpha_background=0.95;   
alpha_defect=1;

% Stacking energies
delta_background=1;     
delta_defect_l=1.5;
delta_defect_r=2;

xi=1E-3;                % Ring factor
c=2.11;                 % loop exponent 
M=1000;                 % number of bps	
i0=499;                 % number of the bp which is defect, should be in the range 1 to M
tolp=1E-9;              % tolerance in root finding method for z_0
z_min=0;                % end-points of the root finding method for z0
z_max=1;

% ----- Notations used in the following ----------

bl=alpha_defect*delta_defect_l;
br=alpha_background*delta_defect_r;
bb=alpha_background*delta_background;
sb=xi*delta_background; 
sl=xi*delta_defect_l;
sr=xi*delta_defect_r;

ab=alpha_background;   
ad=alpha_defect;
db=delta_background;    
dl=delta_defect_l;
dr=delta_defect_r;


% ------------------------------------------------------------
% --------------------- Poland algorithm  --------------------
% -------- Probability that the bp [k,...,k+gammap-1] --------
% --------------- are closed for clamped ends ----------------
% ------------------------------------------------------------  


% --- Assign statistical weights ----------
% --- as determined by the sequence -------


alpha=zeros(M,1);
delta=zeros(M,1);

for i=1:1:M  % loop over bps 
   if i==i0
      alpha(i)=alpha_defect;
      delta(i)=delta_defect_l;
   elseif i==i0+1
       delta(i)=delta_defect_r;
       alpha(i)=alpha_background;
   else
       delta(i)=delta_background;
       alpha(i)=alpha_background;
   end;                   
end;

delta(M+1)=delta(M);

% --- Poland algorithm -------

% assign memory

T=zeros(M+1,1);        
a=zeros(M,1);
a_old=zeros(M,1);
Pu=zeros(M+2,1);
mu=zeros(M+1,1);
mu_old=zeros(M+1,1);


% NOTE: below the arguments in T(k), a(k,j), Pu(k), and mu(k,j) 
% are shifted by one, due to the Matlab indexing convention!

% -- Determine conditional probability [proportional to T_k] --

% Initialize
T(M+1)=alpha(M)*delta(M);

% Iterate
for k=M-1:-1:1
   a(2)=xi*delta(k+2)*g_m(1,c)*T(k+2); 
   for j=2:1:M-k
      a(j+1)=g_m(j,c)*T(k+2)*a_old(j)/g_m(j-1,c);
   end;
   T(k+1)=alpha(k)*delta(k)/(1+sum(a(2:M-k+1)));
   a_old=a;
end;
a(2)=xi*delta(2)*g_m(1,c)*T(2);
for j=2:1:M
   a(j+1)=g_m(j,c)*T(2)*a_old(j)/g_m(j-1,c);
end;
T(1)=1/(1+sum(a(2:M+1)));

% -- Determine unconditional probability  --- 
% -- that bp k (=0,1,...,M+1) is closed   ---

% Initialize
Pu(1)=1;
Pu(2)=T(1);

% Iterate
for k=1:1:M
   if k==1 
      mu(1)=xi*delta(2)*g_m(1,c)*T(1)*T(2);
   else % k>=2delta
      mu(k)=xi*delta(k+1)*g_m(1,c)*T(k)*T(k+1)/(alpha(k-1)*delta(k-1));
   end;
   for j=0:1:k-2
      mu(j+1)=g_m(k-j,c)*T(k+1)*mu_old(j+1)/g_m(k-j-1,c)*delta(k+1)/(delta(k));
   end;
   Pc=T(k+1)/(alpha(k)*delta(k));
   Pu(k+2)=Pu(k+1)*Pc+sum((mu(1:k).*Pu(1:k))); 
   mu_old=mu;
end;

% Check: we must have Pu(M+2)=1
Pu(M+2)



% -------------------------------------------------------------
% ----------------- Analytical solution -----------------------
% -------------------------------------------------------------

z0=bissect_z(z_min,z_max,bb,sb,c,tolp);

Li=(bb-z0)/(sb*z0);            % Li_c(z0), Eq 8
Li_c_1=polylogT(c-1,z0,tolp);  % Li_{c-1}(z0)

P_analy_background=1-1/(1+sb*z0*Li_c_1/bb);  % Eq 9
P_analy_defect=1- 1/(1+bl*sr*Li_c_1/z0*1/(1+sr*Li)*1/(1+sl*Li));  % Eq 6

P=zeros(1,M+2);
D=40;  % proba is calculated in the region [i0-D;i0+D] only
delta=0.0001;  % parameter in the branch cut integral, should be taken close to zero
bc=@(k) branch_cut(delta,k,bb,sb,c);  % bc = I

for i=1:M+1
    i
   if i==i0
       P(i+1)= P_analy_defect;

   elseif (0<(i-i0)) && ((i-i0)<D) % right side of the defect    

         P(i+1)=P_analy_background + z0^(i-i0-1)*(bc(i-i0)*br*(bl-z0*(1+sl*Li))+...   % Eq 12
              bc(i-i0-1)*(sr-sb)/sb*z0*(1+sl*Li))/(bl*sr*Li_c_1+...
              z0*(1+sr*Li)*(1+sl*Li));  

   elseif (0<(i0-i)) && ((i0-i)<D) % left side of the defect  

         P(i+1)=P_analy_background + z0^(i0-i-1)*(bc(i0-i)*br*sl/sr*(bl*sr/sl-z0*(1+sr*Li))+...  
              bc(i0-i-1)*(sl-sb)/sb*z0*(1+sr*Li))/(bl*sr*Li_c_1+...
              z0*(1+sr*Li)*(1+sl*Li));     

    else 
        P(i+1)=P_analy_background; % Further than D, the background proba is used
    end
end



% +++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++ Plots ++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++


d=40;  % We plot the region [i0-d;i0+d] only
k=i0-d:1:i0+d;

figure(1);

% -- Poland algorithm --

plot(k,1-Pu(i0-d:i0+d),'Color',"#80B3FF",'Linewidth',2);  
hold on;

% -- Analytics --

plot(k,P(i0-d:i0+d),'o','MarkerEdgeColor','red','MarkerFaceColor','red','MarkerSize',7);
hold on;

% -- Background proba --

plot([i0-d i0+d],[P_analy_background P_analy_background],'-.k');
hold on;

% -- Defect scope --

d0=-1/log(z0);

plot([i0-d0+1 i0-d0+1],[P_analy_background P_analy_defect],'--k','LineWidth',2);
hold on;
plot([i0+d0+1 i0+d0+1],[P_analy_background P_analy_defect],'--k','LineWidth',2);
hold off;


% ---- Set figure properties ----

L=legend('numerical','analytical','Box','off');
set(L,'Interpreter','latex','Fontsize', 20,'Location', 'northeast');
L.ItemTokenSize(1) = 30;
xlabel('$$k$$','Interpreter','latex','Fontsize', 25);
ylabel('$$1-P(k)$$','Interpreter','latex','Fontsize', 25);

ax=gca; ax.YAxis.Exponent = -3;

set(gca,'DefaultAxesTickLabelInterpreter','latex','Fontsize',15)

xlim([i0-d+2 i0+d])
ylim([min(P(i0-d:i0+d))*0.8 max(P(i0-d:i0+d))*1.1])



