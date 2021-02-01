close all
clear all
a=dlmread ('xx.txt'); % need impedance data in freq. real, imaginary part 
%of impedance with '-' sign
lambda=dlmread('lambda.txt'); % lambda is the regularization parametr
%could be added by user chosen values or from L-curve methods (see below)
f=a(:,1);
z1=a(:,2);
z2=-a(:,3);
Rpol=abs((z1(1)-z1(end)));
h=abs(log10(f(2))-log10(f(1)));
c=f(1:end);
tt=-log10(c*pi*2);
dd=zeros(length(c));
dd2=zeros(length(c));
%z=zeros(length(f));
for i=1:length(c)
    for j=1:length(c)
     %   a2=1/(sqrt(2*pi)*0.33);
        dd(i,j)=(h)/(1+(2*pi*c(i)*10^tt(j))^2)^1;
    %  a3(j)=exp((tt(j)+7.602)/(2*0.33^2));
     % g(j)=a2/a3(j);
       % dd(i,j)=g(j)*(l+1)/(1+(2*pi*f(i)*10^tt(j))^2)^2;
     dd2(i,j)=(h)*(2*pi*c(i)*10^tt(j))/(1+(2*pi*c(i)*10^tt(j))^2)^1;   
    end
end
d=[dd; dd2];  z=[z1;z2];
[U,s,V]=csvd(d);
lambda=l_curve(U,s,z); % the L-curve method could be obtained from 
%"http://www.imm.dtu.dk/~pcha/Regutools/"
 [x_lambda1] = tikhonov(U,s,V,z,lambda); % the Tikhonov regularization 
 %technique could be obtained from "http://www.imm.dtu.dk/~pcha/Regutools/"
 Z1new=dd*x_lambda1;
 Z2new=dd2*x_lambda1;
 x_lambda1=max(x_lambda1,0)/Rpol;
 plot(tt,x_lambda1,'-o');
 xlabel ('log_{10}(\tau)(s)');
 ylabel('DRT');
figure;
plot(z1,z2,'+g');
hold on
 plot(Z1new,Z2new);
 legend('Original','Simulated','Location','Best')
 xlabel ('Z''(\Omega)');
 ylabel('-Z''''(\Omega)');
figure;
ALL=[tt x_lambda1 Z1new Z2new z1 z2];
fileID=fopen('myfile.dat','w');
 save('myfile.dat','ALL','-ascii');
 fclose(fileID);
 