%##################################################
%MODELO DINAMICO USANDO L-E (v2)
%##################################################
close all
clear all
clc

syms L1 L2 m1 m2 g
syms q1 q2 qd1 qd2 qdd1 qdd2

r11=[0;0;0;1]
r22=[0;0;0;1]
g0=[0,0,-g,0]

%% L-E 2.-
A01=[cos(q1),-sin(q1), 0, L1*cos(q1);
     sin(q1), cos(q1), 0, L1*sin(q1);
          0 ,      0 , 1,         0 ;
          0 ,      0 , 0,         1 ];
A12=[cos(q2),-sin(q2), 0, L2*cos(q2);
     sin(q2), cos(q2), 0, L2*sin(q2);
          0 ,       0, 1,         0 ;
          0 ,       0, 0,         1 ];
A02=A01*A12;

%% L-E 3.-
U11=simplify(diff(A01,q1));
U12=simplify(diff(A01,q2));
U21=simplify(diff(A02,q1));
U22=simplify(diff(A02,q2));

%% L-E 4.-
U111=simplify(diff(U11,q1));
U112=simplify(diff(U11,q2));
U121=simplify(diff(U12,q1));
U122=simplify(diff(U12,q2));
U211=simplify(diff(U21,q1));
U212=simplify(diff(U21,q2));
U221=simplify(diff(U22,q1));
U222=simplify(diff(U22,q2));

%% L-E 5.- Obtencion de las matrices de pseudoinercias Ji
%Como se asumio que toda la masa se concentra en el origen
%del sistema de coordenadas, las matrices Ji se reducen a:
J1=[0, 0, 0,  0; 
    0, 0, 0,  0;
    0, 0, 0,  0;
    0, 0, 0, m1];
J2=[0, 0, 0,  0; 
    0, 0, 0,  0;
    0, 0, 0,  0;
    0, 0, 0, m2];

%% L-E 6.- 
%i=1, j=1, k=max(i,j):n => k=1:2
d11=trace(U11*J1*transpose(U11))+trace(U21*J2*transpose(U21));
d11=simplify(d11);

%i=1, j=2, k=max(i,j):n => k=2
d12=trace(U22*J2*transpose(U21));
d12=simplify(d12);

%i=2, j=1, k=max(i,j):n => k=2
d21=trace(U21*J2*transpose(U22));
d21=simplify(d21);

%i=2, j=2, k=max(i,j):n => k=2
d22=trace(U22*J2*transpose(U22));
d22=simplify(d22);

D=[d11,d12;
   d21,d22]

%% L-E 7.- Obtencion las matrices h_ikm
%i=1 k=1 m=1 => j=max(i,k,m):n => j=1:2
h111=trace(U111*J1*transpose(U11))+trace(U211*J2*transpose(U21));
h111=simplify(h111);

%i=1 k=1 m=2 => j=max(i,k,m):n => j=2
h112=trace(U212*J2*transpose(U21));
h112=simplify(h112);

%i=1 k=2 m=1 => j=max(i,k,m):n => j=2
h121=trace(U221*J2*transpose(U21));
h121=simplify(h121);

%i=1 k=2 m=2 => j=max(i,k,m):n => j=2
h122=trace(U222*J2*transpose(U21));
h122=simplify(h122);

%i=2 k=1 m=1 => j=max(i,k,m):n => j=2
h211=trace(U211*J2*transpose(U22));
h211=simplify(h211);

%i=2 k=1 m=2 => j=max(i,k,m):n => j=2
h212=trace(U212*J2*transpose(U22));
h212=simplify(h212);

%i=2 k=2 m=1 => j=max(i,k,m):n => j=2
h221=trace(U221*J2*transpose(U22));
h221=simplify(h221);

%i=2 k=2 m=2 => j=max(i,k,m):n => j=2
h222=trace(U222*J2*transpose(U22));
h222=simplify(h222);

%% L-E 8.- Obtencion de la matriz de fuerzas Centripetas y Coriolis H=[hi]
%i=1 k=1:2 m=1:2
h1=(h111*qd1*qd1+h112*qd1*qd2)+(h121*qd2*qd1+h122*qd2*qd2);
%i=2 k=1:2 m=1:2
h2=(h211*qd1*qd1+h212*qd1*qd2)+(h221*qd2*qd1+h222*qd2*qd2);

H=transpose([h1,h2]);
H=simplify(H)

%% L-E 9.- Obtencion de la matriz de fuerzas de gravedad C=[ci]
%i=1 j=1:2
c1=(-m1*g0*U11*r11)+(-m2*g0*U21*r22);
%i=2 j=1:2
c2=(-m1*g0*U12*r11)+(-m2*g0*U22*r22);

C=transpose([c1,c2]);
C=simplify(C)

%% L-E 10.- La ecuacion dinamica sera:
qdd=[qdd1;qdd2];
tau=D*qdd+H+C
%tau=simplify([tau1;tau2])









