%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Melting curves for Homo-DNA with a defect 
% For background, defect and neighboring basepairs
% For c<1, 1<c<2 and c>2
% Background: A, defect: G
%
% Writes results in txt files, that are then opened to plot
% Comment/uncomment accordingly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
format long;

% import MMT.Core.calculate_breaking_weights;

% ------- Input parameters -------

xi=1E-3;                % Ring factor
saltConc_molar=0.05;    % Molar salt concentration
t_m=65.52;              % melting temparture for AA...AA at saltConc_molar=0.05
tolp=1E-9;              % tolerance in root finding method for z0
z_min=0;                % end-points of the root finding method for z0
z_max=1;
tolT=1E-5;              % tolerance in root finding method for T_star
T_min=t_m+273.15;       % end-points of the root finding method for T_star
T_max=T_min+0.5;   
dis=10;                 % distance to the defect for the neighbour basepair
delta=0.0001;           % used in the branch cut integral, take close to 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------- c<1 -------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=0.5;                 

t_start=63;
t_stop=68; 
t_inc=(t_stop-t_start)/100;

% ------------ Data file -----------------

curve=fopen('melting_curves_05.txt','w');
fclose(curve);

% ------------ Temperature loop -----------------

for t_cel=t_start:t_inc:t_stop
    
    t_cel
    t_kel=t_cel+273.15;   
    [u_hbW_aTaS, u_hbS_aTaS, mat_u_st_ACGT_aTaS] = calculate_breaking_weights(t_kel, saltConc_molar);
    
    ab=u_hbW_aTaS;   
    ad=u_hbS_aTaS;
    db=mat_u_st_ACGT_aTaS(1,1);    
    dl=mat_u_st_ACGT_aTaS(1,3);
    dr=mat_u_st_ACGT_aTaS(3,1);  
    bl=ad*dl;
    br=ab*dr;
    bb=ab*db;
    sb=xi*db; 
    sl=xi*dl;
    sr=xi*dr;
     
    z0=bissect_z(z_min,z_max,bb,sb,c,tolp);
  
    Li=(bb-z0)/(sb*z0);            % Li_c(z0), Eq 8
    Li_c_1=polylogT(c-1,z0,tolp);  % Li_{c-1}(z0)
    
    bc=@(k) branch_cut(delta,k,bb,sb,c);    
    
    % ------------- Probabilities -----------------
    
    P_background=1-1/(1+sb*z0*Li_c_1/bb);  % Eq 9
    
    P_defect=1- 1/(1+bl*sr*Li_c_1/z0*1/(1+sr*Li)*1/(1+sl*Li));  % Eq 6
    
    P_neighbour=P_background + z0^(dis-1)*(bc(dis)*br*(bl-z0*(1+sl*Li))+...   % Eq 12
              bc(dis-1)*(sr-sb)/sb*z0*(1+sl*Li))/(bl*sr*Li_c_1+...
              z0*(1+sr*Li)*(1+sl*Li));  
     
    % --------------- Data export ------------------
    
    curve=fopen('melting_curves_05.txt','a');
         fprintf(curve,'%.5f %.10f %.10f %.10f\n',t_cel,P_background,P_defect,P_neighbour);
    fclose(curve);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------- 1<c<2 -------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=1.5;                

T_star=bissect_T(T_min,T_max,xi,c,saltConc_molar,tolT);
t_star=T_star-273.15; 

t_start=65.4;
t_int=65.5;
t_end=t_star-1e-3;
t_inc_large=(t_int-t_start)/3;
t_inc_small=(t_end-t_int)/40;
   

% --------------- Data file ------------------

curve=fopen('melting_curves_15.txt','w');
fclose(curve);

% ------------ 1st Temperature loop -----------------
% ---------- far from the transition ----------------

for t_cel=t_start:t_inc_large:t_int
    
    t_cel
    t_kel=t_cel+273.15;   
    [u_hbW_aTaS, u_hbS_aTaS, mat_u_st_ACGT_aTaS] = calculate_breaking_weights(t_kel, saltConc_molar);
    
    ab=u_hbW_aTaS;   
    ad=u_hbS_aTaS;
    db=mat_u_st_ACGT_aTaS(1,1);    
    dl=mat_u_st_ACGT_aTaS(1,3);
    dr=mat_u_st_ACGT_aTaS(3,1);  
    bl=ad*dl;
    br=ab*dr;
    bb=ab*db;
    sb=xi*db; 
    sl=xi*dl;
    sr=xi*dr;
    
    z0=bissect_z(z_min,z_max,bb,sb,c,tolp);

    Li=(bb-z0)/(sb*z0);            % Li_c(z0), Eq 8
    Li_c_1=polylogT(c-1,z0,tolp);  % Li_{c-1}(z0)
    
    bc=@(k) branch_cut(delta,k,bb,sb,c);    
    
    % ------------- Probabilities -----------------
    
    P_background=1-1/(1+sb*z0*Li_c_1/bb);  % Eq 9
    
    P_defect=1- 1/(1+bl*sr*Li_c_1/z0*1/(1+sr*Li)*1/(1+sl*Li));  % Eq 6
    
    P_neighbour=P_background + z0^(dis-1)*(bc(dis)*br*(bl-z0*(1+sl*Li))+...   % Eq 12
              bc(dis-1)*(sr-sb)/sb*z0*(1+sl*Li))/(bl*sr*Li_c_1+...
              z0*(1+sr*Li)*(1+sl*Li));  
     
    % --------------- Data export ------------------
    
     curve=fopen('melting_curves_15.txt','a');
          fprintf(curve,'%.5f %.10f %.10f %.10f\n',t_cel,P_background,P_defect,P_neighbour);
     fclose(curve);
    
end

% ------------ 2nd Temperature loop -----------------
% ---------- close to the transition ----------------

for t_cel=t_int:t_inc_small:t_end
    
    t_cel    
    t_kel=t_cel+273.15;    
    [u_hbW_aTaS, u_hbS_aTaS, mat_u_st_ACGT_aTaS] = calculate_breaking_weights(t_kel, saltConc_molar);
    
    ab=u_hbW_aTaS;   
    ad=u_hbS_aTaS;
    db=mat_u_st_ACGT_aTaS(1,1);    
    dl=mat_u_st_ACGT_aTaS(1,3);
    dr=mat_u_st_ACGT_aTaS(3,1);
    bl=ad*dl;
    br=ab*dr;
    bb=ab*db;
    sb=xi*db; 
    sl=xi*dl;
    sr=xi*dr;
    
    z0=bissect_z(z_min,z_max,bb,sb,c,tolp);
   
    Li=(bb-z0)/(sb*z0);            % Li_c(z0), Eq 8
    Li_c_1=polylogT(c-1,z0,tolp);  % Li_{c-1}(z0)
    
    delta=0.000012;       % Numerical errors appear for branch cut near transition, decreased delta for better results
                          % However at delta = 1e-5 also errors for precision of integration

    bc=@(k) branch_cut(delta,k,bb,sb,c);    
    
    % ------------- Probabilities -----------------
    
    P_background=1-1/(1+sb*z0*Li_c_1/bb);  % Eq 9
    
    P_defect=1- 1/(1+bl*sr*Li_c_1/z0*1/(1+sr*Li)*1/(1+sl*Li));  % Eq 6
    
    P_neighbour=P_background + z0^(dis-1)*(bc(dis)*br*(bl-z0*(1+sl*Li))+...   % Eq 12
              bc(dis-1)*(sr-sb)/sb*z0*(1+sl*Li))/(bl*sr*Li_c_1+...
              z0*(1+sr*Li)*(1+sl*Li));  
     
    % --------------- Data export ------------------
    
    curve=fopen('melting_curves_15.txt','a');
         fprintf(curve,'%.5f %.10f %.10f %.10f\n',t_cel,P_background,P_defect,P_neighbour);
    fclose(curve);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------- 2<c -------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=2.5;               

T_star=bissect_T(T_min,T_max,xi,c,saltConc_molar,tolT);
t_star=T_star-273.15;  

t_start=65;
t_end=t_star-1e-5;
t_inc=(t_end-t_start)/60;

% --------------- Data file ------------------
          
 curve=fopen('melting_curves_25.txt','w'); 
 fclose(curve);

% ------------ Temperature loop -----------------

for t_cel=t_start:t_inc:t_end
    
    t_cel
    t_kel=t_cel+273.15;   
    [u_hbW_aTaS, u_hbS_aTaS, mat_u_st_ACGT_aTaS] = calculate_breaking_weights(t_kel, saltConc_molar);
    
    ab=u_hbW_aTaS;   
    ad=u_hbS_aTaS;
    db=mat_u_st_ACGT_aTaS(1,1);    
    dl=mat_u_st_ACGT_aTaS(1,3);
    dr=mat_u_st_ACGT_aTaS(3,1);  
    bl=ad*dl;
    br=ab*dr;
    bb=ab*db;
    sb=xi*db; 
    sl=xi*dl;
    sr=xi*dr;
    
    z0=bissect_z(z_min,z_max,bb,sb,c,tolp);

    Li=(bb-z0)/(sb*z0);            % Li_c(z0), Eq 8
    Li_c_1=polylogT(c-1,z0,tolp);  % Li_{c-1}(z0)
    
    delta=0.0001;
    bc=@(k) branch_cut(delta,k,bb,sb,c);    
    
    % ------------- Probabilities -----------------
    
    P_background=1-1/(1+sb*z0*Li_c_1/bb);  % Eq 9
    
    P_defect=1- 1/(1+bl*sr*Li_c_1/z0*1/(1+sr*Li)*1/(1+sl*Li));  % Eq 6
    
    P_neighbour=P_background + z0^(dis-1)*(bc(dis)*br*(bl-z0*(1+sl*Li))+...   % Eq 12
              bc(dis-1)*(sr-sb)/sb*z0*(1+sl*Li))/(bl*sr*Li_c_1+...
              z0*(1+sr*Li)*(1+sl*Li)); 

    % --------------- Data export ------------------
    
     curve=fopen('melting_curves_25.txt','a');
          fprintf(curve,'%.5f %.10f %.10f %.10f\n',t_cel,P_background,P_defect,P_neighbour);
     fclose(curve);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------ Plots ---------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------- c<1 -------

subplot(1,3,1)

data_cl1=importdata('melting_curves_05.txt');

plot(data_cl1(:,1), data_cl1(:,2), data_cl1(:,1), data_cl1(:,3), data_cl1(:,1), data_cl1(:,4),'Linewidth',2)
hold on;

xlabel('$T(^{\circ}C)$','Interpreter','latex','Fontsize', 20);
ylabel('$1-P(k)$','Interpreter','latex','Fontsize', 20);
set(gca,'Fontsize',15)
title('$(a) \ c=0.5$','Interpreter','latex','Fontsize', 20)


% ------- 1<c<2 -------

subplot(1,3,2)

data_1lcl2=importdata('melting_curves_15.txt');

plot(data_1lcl2(:,1), data_1lcl2(:,2),data_1lcl2(:,1), data_1lcl2(:,3),data_1lcl2(:,1), data_1lcl2(:,4),'Linewidth',2)
hold on;

xlabel('$T(^{\circ}C)$','Interpreter','latex','Fontsize', 20);
set(gca,'Fontsize',15)
title('$(b) \ c=1.5$','Interpreter','latex','Fontsize', 20)

M=max([max(data_1lcl2(:,2)); max(data_1lcl2(:,3)); max(data_1lcl2(:,4))]);
m=min([min(data_1lcl2(:,2)); min(data_1lcl2(:,3)); min(data_1lcl2(:,4))]);
axis([65.4 65.58 m*0.8 M*1.2])

% T_star asymptot

c=1.5;                
T_star=bissect_T(T_min,T_max,xi,c,saltConc_molar,tolT);
t_star=T_star-273.15; 
plot([t_star t_star],[0 1],'--k');
hold on;


% ------- 2<c -------

subplot(1,3,3)

data_2lc=importdata('melting_curves_25.txt');

plot(data_2lc(:,1), data_2lc(:,2),data_2lc(:,1), data_2lc(:,3),data_2lc(:,1), data_2lc(:,4),'Linewidth',2);
hold on;

xlabel('$T(^{\circ}C)$','Interpreter','latex','Fontsize', 20);
set(gca,'Fontsize',15)
title('$(c) \ c=2.5$','Interpreter','latex','Fontsize', 20)

M=max([max(data_2lc(:,2)); max(data_2lc(:,3)); max(data_2lc(:,4))]);
m=min([min(data_2lc(:,2)); min(data_2lc(:,3)); min(data_2lc(:,4))]);
axis([65 65.65 m*0.8 M*1.2])

% T_star asymptot

c=2.5;                
T_star=bissect_T(T_min,T_max,xi,c,saltConc_molar,tolT);
t_star=T_star-273.15; 
plot([t_star t_star],[0 1],'--k');
hold on;

% Values B

T_star_eps=1e-5;
[u_hbW_aTaS, u_hbS_aTaS, mat_u_st_ACGT_aTaS] = calculate_breaking_weights(T_star-T_star_eps, saltConc_molar);
bb=u_hbW_aTaS*mat_u_st_ACGT_aTaS(1,1);
bl=u_hbS_aTaS*mat_u_st_ACGT_aTaS(1,3);
br=u_hbW_aTaS*mat_u_st_ACGT_aTaS(3,1);
sb=xi*mat_u_st_ACGT_aTaS(1,1);
sl=xi*mat_u_st_ACGT_aTaS(1,3);
sr=xi*mat_u_st_ACGT_aTaS(3,1);  

delta=0.0001;
bc=@(k) branch_cut(delta,k,bb,sb,c);
   
B_bg=1-1/(1+sb*zeta(c-1)/bb);
B_df=1-1/(1+bl*sr*zeta(c-1)/((1+sr*zeta(c-1))*(1+sl*zeta(c-1))));
B_nb=B_bg+(bc(dis)*br*sb*(bl-(1+sl*zeta(c)))+bc(dis-1)*(sr-sb)*(1+sl*zeta(c)))/(sb*bl*sr*zeta(c-1)+sb*(1+sl*zeta(c))*(1+sr*zeta(c)));

plot(t_star, B_bg,'o','MarkerEdgeColor','#0072BD','MarkerFaceColor','#0072BD','MarkerSize',7);
plot(t_star, B_df,'o','MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319','MarkerSize',7);
plot(t_star, B_nb,'o','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120','MarkerSize',7);


% ------- Subplots and legend positions -------

ha=get(gcf,'children');
set(ha(1),'position',[.682 .12 .266 0.8])
set(ha(2),'position',[.366 .12 .266 0.8])
set(ha(3),'position',[.05 .12 .266 0.8])

L=legend('$k \to \infty$','$k=k_d$','$k=k_d+10$','','Box','off');
set(L,'Interpreter','latex','Fontsize', 20,'Position',[0.3 0.64 0.3 0.25]);

x0=10;
y0=10;
width=2000;
height=350;
set(gcf,'position',[x0,y0,width,height])
set(gcf, 'PaperSize', [40 40])

