clear;
clc;

T = 300;
m0=9.11*10^-31; %in kg
mn=0.26*m0;
kB=1.38*10^-23;
vth=sqrt(2*kB*T/mn);%thermal velocity
deltaV=0.1;
W=1*10^-7;
L=2*10^-7;
E=deltaV/L; %Question 1 (a)
q=1.6*10^-19;
F=q*E; %Question 1 (b)
a=F/mn; %Question 1 (c)
N=20; %number of particles
deltaT=5*10^-15;
tTotal=1000*deltaT;
tmn=0.2*10^-12;
n=10^15;%electron concentration 

%Initial velocity 
vx=zeros(1,N);
vy=zeros(1,N);

%Initial position of each particles
xPos=zeros(1,N);
yPos=zeros(1,N);

for i = 1:N
    x=rand*L;
    y=rand*W;
    xPos(i)=xPos(i)+x;
    yPos(i)=yPos(i)+y;
    vx(i) = vth/sqrt(2)*randn;
    vy(i) = vth/sqrt(2)*randn; 
end


for t = 0 : deltaT : tTotal 
    for i=1:N
        P = 1-exp(-deltaT/tmn);
        if P > rand()
            vx(i) = vth/sqrt(2)*randn;
            vy(i) = vth/sqrt(2)*randn;
        end
    end
   
    deltaPx = vx*deltaT+(1/2)*a*deltaT^2;
    deltaPy = vy*deltaT;
    for i=1:N
        if yPos(i)+deltaPy(i)>W||yPos(i)+deltaPy(i)<0
            vy(i)=-vy(i);
            deltaPy(i)=vy(i)*deltaT;
        end
        
    end
    
    xPos=xPos+deltaPx;
    yPos=yPos+deltaPy;
    vxAvg=mean(vx);
    j=n*q*vxAvg;
    I=j*W
    figure(1);
    plot(t,I,'.');
    xlim([0 tTotal]);
    xlabel('t');
    ylabel('I');
    hold on;
    
    %Periodic boundary condition in x direction
    Ix=xPos>L;
    xPos(Ix)=xPos(Ix)-L;
    Ix=xPos<0;
    xPos(Ix)=xPos(Ix)+L;
    
    figure(2);
    plot(xPos,yPos,'.');
    hold on;
    xlim([0 L]);
    ylim([0 W]);

    KEsum=0;
    for i = 1:N
         v_Squared = vx(i)^2+vy(i)^2;
         KEsum = KEsum + (1/2)*mn*v_Squared;
    end
    KEavg = KEsum /N; 
    T=KEavg/kB;
    
    pause(0.001);
    vx=vx+a*deltaT;

end

%Part3 question c: electron density map 
P=zeros(200,100);
xPos = xPos.*10^9;
yPos = yPos.*10^9;
for i=1:N
    for j = 1:200
        for k = 1:100
            if xPos(i) > j && xPos(i)< (j+1) && yPos(i)>k && yPos(i)<(k+1)
                P(j,k) = P(j,k)+1;
            end
        end
    end
end
figure(3);
surf(P);
title('Electron density map');

%Part3 question d: temperature map 
Temp=zeros(200,100);
xPos = xPos;
yPos = yPos;
for i=1:N
    for j = 1:200
        for k = 1:100
            if xPos(i) > j && xPos(i)< (j+1) && yPos(i)>k && yPos(i)<(k+1)
                v_Squared = vx(i)^2+vy(i)^2;
                T = (1/2)*mn*v_Squared/kB;
                Temp(j,k) = Temp(j,k)+T; 
            end
        end
    end
    KE=0;
end
figure(4); 
surf(Temp);
title('Temperature map');          
hold off;


