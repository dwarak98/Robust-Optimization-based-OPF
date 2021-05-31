

%
clear;
clc;

op=runpf('case30');

%op.branch(1,PT);
% generator data
         % bus no. Energy  pmax    pmin   Rampup   Rampd    AGC    Amax Amin Reserve 
         %       (yuan/MW)  MW     MW    MW/5min  MW/5min (yuan/MW) MW   MW  (yuan/MW)
gen_data=[
            1       150     80      0       12      -12     120     75   15   120;
            2       180     80      0       16      -16     200     70   10   150;
            22      200     50      0       15      -15     180     50   10   100;
            23      250     30      0       8       -8      0       0    0    80;
            27      300     55      0       10      -10     0       0    0    90];

mpc=case_ieee30_0();
H=makePTDF(mpc);
Sa=mpc.branch(:,6);
Sb=mpc.branch(:,7);
Sc=mpc.branch(:,8);
for i=1:length(Sa)
    Fmax(i)=max([Sa(i),Sb(i),Sc(i)]);
end

Fmax=transpose(Fmax);
Gp=[H(:,1) H(:,2) H(:,22) H(:,23) H(:,27)];
Gd=H;
Gw=[H(:,13) H(:,20)];
l=length(Sa);
D=mpc.bus(:,3);
E(1,:)=[-4 4]*0;
Emax(1,1)=max(E(1,:));
Emin(1,1)=min(E(1,:));
E(2,:)=[-8 8]*0;
Emax(2,1)=max(E(2,:));
Emin(2,1)=min(E(2,:));

demand=sum(mpc.bus(:,3));
%W=zeros(30,1);
W(1,1)=20;
W(2,1)=40;
D_spin=30;
D_A=20;
eps=0.8;
cvx_begin

    variables P(5) A_pos(3) A_neg(3) spin(5) BET(3) z1(3) z2(3) y_pos(l,2) y_neg(l,2)
    %variable u(5) binary
    dual variables mu1 mu2
    %P=P.*u;

    minimize (sum(P.*gen_data(:,2)) + sum((A_pos+A_neg).*gen_data(1:3,7)) + sum(spin.*gen_data(:,10)))
    
    subject to
    
        %-------------Load balancing-----------------%
        
        sum(P)==demand-sum(W)
        
        %--------spinning reserve cpnstraint--------%
        
        sum(spin)>=D_spin
        spin >= 0
        
        %---------regulating service ---------------%
        
        sum(A_pos)>=D_A
        sum(A_neg)>=D_A
        P(1:3)+A_pos<= gen_data(1:3,8)
        P(1:3)-A_neg>= gen_data(1:3,9)
        0<=A_pos<=gen_data(1:3,5)
        0<=A_neg<=abs(gen_data(1:3,6))
        %A_pos==A_neg
        % ramp up constraint for spin reserve
        %--------- Generation dispatch limits--------%
        
        P(1:3)+A_pos(1:3)+spin(1:3)<=gen_data(1:3,3)
        P(1:3)-A_neg(1:3)+spin(1:3)>=gen_data(1:3,4)
        
        gen_data(4:5,4)<=P(4:5)+spin(4:5)<=gen_data(4:5,3)
        P>=0
        
        %-----------transmission line constraints------%
        %
        y_pos>=0
        y_pos>= [Gp(:,1:3)*BET,Gp(:,1:3)*BET] + Gw
        
        y_neg<=0
        y_neg<= [Gp(:,1:3)*BET,Gp(:,1:3)*BET] + Gw
        
        
        
        
        
        mu1: Gp*P-Gd*D+Gw*W+y_pos*Emax+y_neg*Emin <= Fmax*eps 
        
        mu2: -Fmax*eps <= Gp*P-Gd*D+Gw*W+y_pos*Emin+y_neg*Emax
        
        %-------------uncertainty constraints----------%
        sum(BET)==1;
        %BET>=0
        
        z1>=0
        z1>=BET
        z2<=0
        z2<=BET
        z1*sum(Emax)+z2*sum(Emin)<=A_pos
        z1*sum(Emin)+z2*sum(Emax)>=-A_neg
        %}
        
        
        
 
   
cvx_end
A_pos(4:5)=0;
A_neg(4:5)=0;
bus_no=gen_data(:,1);
T=table(bus_no,P,A_pos,A_neg,spin)

        





    