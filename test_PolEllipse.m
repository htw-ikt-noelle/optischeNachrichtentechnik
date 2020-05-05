% 2020-04-12
% see Master-Thesis Mannheimer2014, p.6 / Fig.2.2

w=1.5*pi*1; 
T=2.0; 
Fs=200*1/T; 
t=0:1/Fs:T-1/Fs;
A = 1.4;

Ax=1; phix=-pi/2*1; 
X=Ax*cos(w*t+phix); 

Ay=1; phiy=0; 
Y=Ay*cos(w*t+phiy); 

z=linspace(0,1,length(t));

figure(1),clf
plot3(t,X,z*0-A,'--', t,z*0+A,Y,'--', t,X,Y, z*0,X,Y, t,z*0,z*0,'k','LineWidth',1.5), 
ax = gca;

grid on, 
axis([0 T -1 1 -1 1]*A), axis square, 
xlabel('x'),ylabel('y'),zlabel('z')
view(45,20)
