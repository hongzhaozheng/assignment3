clear;
clc;


W=100;
L=200;
V0=0.1;
so=1;
si=10^-2;
Wb=(2/5)*W;
Lb=(2/5)*L;

%block definition
line([Lb Lb],[0 -Wb]);
line([Lb Lb],[-W Wb-W]);
line([2*Lb L-Lb],[0 -Wb]);
line([2*Lb L-Lb],[-W Wb-W]);
line([Lb L-Lb],[-Wb -Wb]);
line([Lb L-Lb],[Wb-W Wb-W]);


for i=1:W
    for j=1:L
        n = i+(j-1)*W;
        if (j>Lb && j<L-Lb) && (i<Wb || i>(W-Wb))
            sigma(i,j)=si;
        else
            sigma(i,j)=so;
        end
        sigmaV(n) = sigma(i,j);
    end
end

G=sparse(W*L,W*L);

%n:i,j
%nl:i,j-1
%nr:i,j+1
for i=1:(W*L)
    if i <= W || i > (L-1)*W
        G(i,i)=1;
    end
end

for j=2:(L-1)
    for i=1:W
        n = i+(j-1)*W;
        a= i+(j-2)*W;
        d= i+j*W;
        w= (i-1)+(j-1)*W;
        s= (i+1)+(j-1)*W;
        
        if i==1
            G(n,s)=1/(1/sigma(n) + 1/sigma(s));
            G(n,a)=1/(1/sigma(n) + 1/sigma(a));
            G(n,d)=1/(1/sigma(n) + 1/sigma(d));
            G(n,n)=-(G(n,s) + G(n,a) + G(n,d));
        elseif i==W
            G(n,w)=1/(1/sigma(n) + 1/sigma(w));
            G(n,a)=1/(1/sigma(n) + 1/sigma(a));
            G(n,d)=1/(1/sigma(n) + 1/sigma(d));
            G(n,n)=-(G(n,w) + G(n,a) + G(n,d));
        else
            G(n,s)=1/(1/sigma(n) + 1/sigma(s));
            G(n,w)=1/(1/sigma(n) + 1/sigma(w));
            G(n,a)=1/(1/sigma(n) + 1/sigma(a));
            G(n,d)=1/(1/sigma(n) + 1/sigma(d));
            G(n,n)=-(G(n,w) + G(n,a) + G(n,d) + G(n,s));
        end
        
    end
end

F=sparse(W*L,1);

for i=1:(W*L)
    if i <= W
        F(i,1)=V0;
    end
end

V=G\F;

for i=1:W
    for j=1:L
        n=i+(j-1)*W;
        Emap(i,j) = V(n);
    end
end
figure(1);
mesh(Emap);
title('V(x,y)');

[X Y]=gradient(Emap);

Ex=-X;
Ey=-Y;

figure(2);
quiver(Ex,Ey);
title('E(x,y)');
