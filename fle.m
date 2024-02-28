%##############################################################
% MODELO DINAMICO L-E(q,qd,qdd)
%Este codigo halla el modelo dinamico de un robot de 2 GDL
%que trabaja en el plano X,Y.
%El modelo devuelto es una expresion simbolica/Ecuacion literal
%##############################################################
close all
clear all
clc

%parametros cinematicos y dinamicos
syms m1 m2 L1 L2 g

%variables articulares
syms q1 q2
syms qd1 qd2
syms qdd1 qdd2

%% Condiciones:
%-Asumimos que toda la masa del eslabon i se concentra en
% el origen del sistema {Si} asociado, es decir, en el extremo
% final del eslabon.
r11=[0;0;0;1]
r22=[0;0;0;1]
%-Asumimos la direccion de la gravedad en z y en sentido negativo
% respecto al sistema de la base {S0}.
g0=[0,0,-g,0]

%% L-E 2.- MTH intermedias
%Previamente se hallo la tabla de parametros D-H:
%-------------------------------
%|  tethai | di  | ai | alphai |
%|------------------------------
%|    q1   |  0  | L1 |   0    |
%|    q2   |  0  | L2 |   0    |
%-------------------------------

A01=[cos(q1),-sin(q1), 0, L1*cos(q1);
     sin(q1), cos(q1), 0, L1*sin(q1);
          0 ,      0 , 1,         0 ;
          0 ,      0 , 0,         1 ];
A12=[cos(q2),-sin(q2), 0, L2*cos(q2);
     sin(q2), cos(q2), 0, L2*sin(q2);
          0 ,       0, 1,         0 ;
          0 ,       0, 0,         1 ];

A02=A01*A12;

%% L-E 3.- Obtencion de las matrices Uij
%Uij=d(Ai)/dqj
%En lugar de derivar parcialmente, se utiliza la expresion
% sugerida en la nota 1 en el libro.
Q1=[0, -1, 0, 0;
    1,  0, 0, 0;
    0,  0, 0, 0;
    0,  0, 0, 0];
Q2=[0, -1, 0, 0;
    1,  0, 0, 0;
    0,  0, 0, 0;
    0,  0, 0, 0];
M0=[0, 0, 0, 0;   %matriz nula
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0];

%Calculo de las matrices Uij:
U11=Q1*A01;      %i=1  j=1
U12=M0;          %i=1  j=2
U21=Q1*A01;      %i=2  j=1
U22=A01*Q2*A12;  %i=2  j=2

%% L-E 4.- Obtencion de las matrices Uijk
%Uijk=d(Uij)/dqk
%En lugar de derivar parcialmente, se utiliza la expresion
% sugerida en la nota 2 en el libro.

U111=Q1*Q1*A01;      %i=1  j=1  k=1
U112=M0;             %i=1  j=1  k=2
U121=M0;             %i=1  j=2  k=1
U122=M0;             %i=1  j=2  k=2
U211=Q1*Q1*A02;      %i=2  j=1  k=1
U212=Q1*A01*Q2*A12;  %i=2  j=1  k=2
U221=Q1*A01*Q2*A12;  %i=2  j=2  k=1
U222=A01*Q2*Q2*A12;  %i=2  j=2  k=2

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

%% L-E 6.- Obtencion de la matriz de inercias D=[dij]

%i=1, j=1, k=max(i,j):n => k=1:2
d11=traza(U11*J1*transpose(U11))+traza(U21*J2*transpose(U21));
d11=simplify(d11);

%i=1, j=2, k=max(i,j):n => k=2
d12=traza(U22*J2*transpose(U21));
d12=simplify(d12);

%i=2, j=1, k=max(i,j):n => k=2
d21=traza(U21*J2*transpose(U22));
d21=simplify(d21);

%i=2, j=2, k=max(i,j):n => k=2
d22=traza(U22*J2*transpose(U22));
d22=simplify(d22);

D=[d11,d12;
   d21,d22]

%% L-E 7.- Obtencion las matrices h_ikm
%i=1 k=1 m=1 => j=max(i,k,m):n => j=1:2
h111=traza(U111*J1*transpose(U11))+traza(U211*J2*transpose(U21));
h111=simplify(h111);

%i=1 k=1 m=2 => j=max(i,k,m):n => j=2
h112=traza(U212*J2*transpose(U21));
h112=simplify(h112);

%i=1 k=2 m=1 => j=max(i,k,m):n => j=2
h121=traza(U221*J2*transpose(U21));
h121=simplify(h121);

%i=1 k=2 m=2 => j=max(i,k,m):n => j=2
h122=traza(U222*J2*transpose(U21));
h122=simplify(h122);

%i=2 k=1 m=1 => j=max(i,k,m):n => j=2
h211=traza(U211*J2*transpose(U22));
h211=simplify(h211);

%i=2 k=1 m=2 => j=max(i,k,m):n => j=2
h212=traza(U212*J2*transpose(U22));
h212=simplify(h212);

%i=2 k=2 m=1 => j=max(i,k,m):n => j=2
h221=traza(U221*J2*transpose(U22));
h221=simplify(h221);

%i=2 k=2 m=2 => j=max(i,k,m):n => j=2
h222=traza(U222*J2*transpose(U22));
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

%% Funcion para hallar la traza de una matriz 4x4
function res = traza(X)
%Traza: es la suma de los elementos de la diagonal principal de una matriz
res=X(1,1)+X(2,2)+X(3,3)+X(4,4);
end

