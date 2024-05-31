clear
close all;
clc

syms al be ep u real
assume((0<=ep)&(ep<=1))
assume((0<=al)&(al<=1))
assume((0<=be)&(be<=1))
assume((0<=u)&(u<=1))
assume(al^2+be^2==1)
assume(al,'real')
assume(be,'real')
assume(u,'real')

%Initial state parameters:
H = [1;0]; %|0>
V = [0;1]; %|1>
zero=H;
Id2 = eye(2,2);

%Purification circuit parameter and operations:
aqubit=[1,0;0,0];
CNOT=[1,0,0,0;0,1,0,0;0,0,0,1;0,0,1,0];
CNOT_R=[1,0,0,0;0,0,0,1;0,0,1,0;0,1,0,0];

kk=0;
for b=0:0.005:1
    a=sqrt(1-b^2);
    kk=kk+1;
    jj=0;
    for la=0:0.005:1
        jj=jj+1;
        si=sqrt(1-b^2)*kron(H,H)+b*kron(V,V);
        rnd=(1-la)*(1/4)*kron(Id2,Id2)+(la)*(si*si');
        rndn=rnd/trace(rnd);
        conc(kk,jj)=  max(0,(2*(abs(rndn(1,4) )-sqrt(rndn(2,2)*rndn(3,3)))));       
        Con_cal=max (0,2*la*a*b - abs(la - 1)/2);
        th=atan((sqrt(1-b^2)/b));
        H1=[cos(th),-sin(th);sin(th),cos(th)];
        rho_a1=kron(CNOT_R,Id2)*kron(H1*aqubit*H1',rndn)*kron(CNOT_R,Id2)';
        M0=kron(zero*zero',kron(Id2,Id2));
        M1=kron(V*V',kron(Id2,Id2));
        rho_M0=M0*rho_a1*M0';
        rhof=(PartialTrace(rho_M0,[1]));
        s0=trace(rhof);
        rhofn0=(rhof/s0);
        conm(kk,jj)=max(0, (2*(abs(rhofn0(1,4) )-sqrt(rhofn0(2,2)*rhofn0(3,3)))));
        %Difference of concurrence:
        cdiff=conm(kk,jj)-conc(kk,jj);
        if cdiff<0
            cd1(kk,jj)=-1;
            cd2(kk,jj)=NaN;
            cd3(kk,jj)=NaN;
        else if cdiff==0
               cd2(kk,jj)=0;
               cd1(kk,jj)=NaN;
               cd3(kk,jj)=NaN;
            else cdiff>0
               cd3(kk,jj)=1;
               cd1(kk,jj)=NaN;
               cd2(kk,jj)=NaN;
            end
        end
       
    end
end
la=0:0.005:1;
et=0:0.005:1;
figure(1)
surf(et,la,conc) 
axis tight
xlabel('F')
ylabel('\beta')
zlabel('C(\rho^{dp})')
colormap(hot)
colorbar
 
figure(2)
surf(et,la,conm) 
axis tight
xlabel('F')
ylabel('\beta')
zlabel('C(\rho^{dp}_{PC(opt)})')
colormap(hot)
colorbar
    
figure(3)
surf(cd1); hold on
surf(cd2);hold on
surf(cd3)
axis tight
xlabel('F')
ylabel('\beta')
box on;
