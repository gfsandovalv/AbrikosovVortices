%------------------------------------------------------------------------%
%--SIMULACION MODELO LINEAL G-L INDEPENDIENTE DEL TIEMPO-----------------%
%---------- Juan Pablo Herrera, Eliana Vargas, Camilo Fajardo, Gabriel Sandoval, Santiago Meléndez-------------%
%------------------------------------------------------------------------%
clc
clf
close all
clear N D2 V2 Ki Kd C C0 C1 C2 M1 A d Mo1 Mo K f
clear all
%--------------------------------------------------------------------------
%-------------------- INGRESO DE DATOS INICIALES------------
%--------------------------------------------------------------------------
so=1;
co=10000000000;
Ks=180.3;                                %Constante K modelo lineal
He=0.000003;                          %Campo externo aplicado
%-------------------------------------------------------------------------%
%--------------------LECTURA ARCHIVO DE MALLADO GMSH-----------------------
%-------------------------------------------------------------------------%

file = ('cell.msh');%-----------archivo txt creado por gmsh------------------
N_n= dlmread(file,'',[5-1 1-1 5-1 1-1]);
N_e= dlmread(file,'',[7+N_n 0 7+N_n 0]);
node_id = dlmread(file,'',[5 0 4+N_n 0]);
nodes = dlmread(file,'',[5 1 4+N_n 3]);
elements = dlmread(file,'',[8+N_n 0 7+N_n+N_e 7]);
%------- Geometría 2D-----------------------------------------------------
two_d_nodes = nodes(:,1:2);
elem_type = elements(:,2);
%--- Halla el índice inicial de los elementos 2D----------------------
two_ind = 1;
for i = 1:N_e
    if(elem_type(i) ~= 2)
        two_ind = two_ind+1;
    end
end
%----------------------------------------------
two_d_elements(1:N_e-two_ind,1:3) = 0;
k = 1;
for i = two_ind:N_e
    two_d_elements(k,1:3) = elements(i,6:8);
    k = k+1;
end
%---------Asignacion de código a elementos de frontera------------

one_d_elements=two_ind-1;
Bound=zeros(one_d_elements,3);
for i=1:1:one_d_elements
    Bound(i,1)=elements(i,4);
    Bound(i,2)=elements(i,6);
    Bound(i,3)=elements(i,7);
end
%--------------------Cálculo de pendientes en frontera---------------

N_bound=0;
n=1;
s=0;
i=1;
while i~=one_d_elements+1
    j=n;
    while elements(i,4)-elements(j,4)==0
        n=n+1;
        j=j+1;
    end
    i=n;
    N_bound=N_bound+1;  
end
%------Calcula número de fronteras
Pend=zeros(N_bound,1);
Codigo=zeros(N_bound,1);
Vec=zeros(N_bound,2);
k=1;
i=1;
n=1;
while i~=one_d_elements+1
    j=n;
    while elements(i,4)-elements(j,4)==0
        n=n+1;
        j=j+1;
    end
    Codigo(k,1)=elements(j-1,4);
    Pend(k,1)=(nodes(Bound(i,3),2)-nodes(Bound(i,2),2))/(nodes(Bound(i,3),1)-nodes(Bound(i,2),1));
    Vec(k,1)=nodes(Bound(i,3),1)-nodes(Bound(i,2),1); %---Normal comp x
    Vec(k,2)=nodes(Bound(i,3),2)-nodes(Bound(i,2),2); %---Normal comp y
    k=k+1;
    i=n;
end
%---------Vector Normal Unitario en Frontera-----------------------

Nor_u=zeros(N_bound,2);

for i=1:1:N_bound
    if Vec(i,2)==0
        Nor_u(i,1)=0;
        Nor_u(i,2)=1;
    else
        Nor_u(i,1)=1;
        Nor_u(i,2)=-Vec(i,1)/Vec(i,2);
    end
    Nor_u(i,:)=Nor_u(i,:)/(Nor_u(i,1).^2+Nor_u(i,2).^2).^0.5;
end

Matriz=zeros(N_bound,4);
Matriz(:,1)=Codigo;
Matriz(:,2)=Pend;
Matriz(:,3)=Nor_u(:,1);
Matriz(:,4)=Nor_u(:,2);

%---- Gráfica de Mallado-------------------------------------------------
figure(1)
triplot(two_d_elements(:,1:3),two_d_nodes(:,1),two_d_nodes(:,2),'k')
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
title('Discretización del Dominio','fontsize',14)
fh = figure(1);
set(fh, 'color', 'white');
hold on
x=nodes(34,1);
y=nodes(34,2);
plot(x,y,'*');
%--------------------------------------------------------------------------
%----------DETERMINACION VECTOR POTENCIAL MAGNÉTICO------------------------
%--------------------------------------------------------------------------
d = ones(3,3);
C0=zeros(6,1);
C1=zeros(6,1);
B=zeros(6,6);
Div=zeros(6,6);
K=sparse(2*N_n,2*N_n);
F=sparse(2*N_n,1);
tic
for i = 1:1:N_e-one_d_elements;
    n1 = two_d_elements(i,1);
    n2 = two_d_elements(i,2);
    n3 = two_d_elements(i,3);
    d(1,2) = nodes(n1,1); d(1,3) = nodes(n1,2);
    d(2,2) = nodes(n2,1); d(2,3) = nodes(n2,2);
    d(3,2) = nodes(n3,1); d(3,3) = nodes(n3,2);
    J =det(d);

 %-----Jacobiano por elemento-------------------
%-------------------Cálculo de Matriz B y Div--------------------------
    C0(1,1)=-(nodes(n3,1)-nodes(n2,1))/J;
    C0(2,1)=-(nodes(n1,1)-nodes(n3,1))/J;
    C0(3,1)=-(nodes(n2,1)-nodes(n1,1))/J;
    C0(4,1)=(nodes(n2,2)-nodes(n3,2))/J;
    C0(5,1)=(nodes(n3,2)-nodes(n1,2))/J;
    C0(6,1)=(nodes(n1,2)-nodes(n2,2))/J;
    B=2*C0*transpose(C0)*J*0.5;
    f=2*He*C0*J*0.5;
    C1(1,1)=(nodes(n2,2)-nodes(n3,2))/J;
    C1(2,1)=(nodes(n3,2)-nodes(n1,2))/J;
    C1(3,1)=(nodes(n1,2)-nodes(n2,2))/J;
    C1(4,1)=(nodes(n3,1)-nodes(n2,1))/J;
    C1(5,1)=(nodes(n1,1)-nodes(n3,1))/J;
    C1(6,1)=(nodes(n2,1)-nodes(n1,1))/J;
    Div=2*so*C1*transpose(C1)*J*0.5;
    %-------------------Ensamblaje Matriz de Rigidez--------------------
    A=B+Div;
    K(n1,n1)=K(n1,n1)+A(1,1);
    K(n1,n2)=K(n1,n2)+A(1,2);
    K(n1,n3)=K(n1,n3)+A(1,3);
    K(n1,n1+N_n)=K(n1,n1+N_n)+A(1,4);
    K(n1,n2+N_n)=K(n1,n2+N_n)+A(1,5);
    K(n1,n3+N_n)=K(n1,n3+N_n)+A(1,6);
    K(n2,n1)=K(n2,n1)+A(2,1);
    K(n2,n2)=K(n2,n2)+A(2,2);
    K(n2,n3)=K(n2,n3)+A(2,3);
    K(n2,n1+N_n)=K(n2,n1+N_n)+A(2,4);
    K(n2,n2+N_n)=K(n2,n2+N_n)+A(2,5);
    K(n2,n3+N_n)=K(n2,n3+N_n)+A(2,6);
    K(n3,n1)=K(n3,n1)+A(3,1);
    K(n3,n2)=K(n3,n2)+A(3,2);
    K(n3,n3)=K(n3,n3)+A(3,3);
    K(n3,n1+N_n)=K(n3,n1+N_n)+A(3,4);
    K(n3,n2+N_n)=K(n3,n2+N_n)+A(3,5);
    K(n3,n3+N_n)=K(n3,n3+N_n)+A(3,6);
    K(n1+N_n,n1)=K(n1+N_n,n1)+A(4,1);
    K(n1+N_n,n2)=K(n1+N_n,n2)+A(4,2);
    K(n1+N_n,n3)=K(n1+N_n,n3)+A(4,3);
    K(n1+N_n,n1+N_n)=K(n1+N_n,n1+N_n)+A(4,4);
    K(n1+N_n,n2+N_n)=K(n1+N_n,n2+N_n)+A(4,5);
    K(n1+N_n,n3+N_n)=K(n1+N_n,n3+N_n)+A(4,6);
    K(n2+N_n,n1)=K(n2+N_n,n1)+A(5,1);
    K(n2+N_n,n2)=K(n2+N_n,n2)+A(5,2);
    K(n2+N_n,n3)=K(n2+N_n,n3)+A(5,3);
    K(n2+N_n,n1+N_n)=K(n2+N_n,n1+N_n)+A(5,4);
    K(n2+N_n,n2+N_n)=K(n2+N_n,n2+N_n)+A(5,5);
    K(n2+N_n,n3+N_n)=K(n2+N_n,n3+N_n)+A(5,6);
    K(n3+N_n,n1)=K(n3+N_n,n1)+A(6,1);
    K(n3+N_n,n2)=K(n3+N_n,n2)+A(6,2);
    K(n3+N_n,n3)=K(n3+N_n,n3)+A(6,3);
    K(n3+N_n,n1+N_n)=K(n3+N_n,n1+N_n)+A(6,4);
    K(n3+N_n,n2+N_n)=K(n3+N_n,n2+N_n)+A(6,5);
    K(n3+N_n,n3+N_n)=K(n3+N_n,n3+N_n)+A(6,6);
    F(n1,1)=F(n1,1)+f(1,1);
    F(n2,1)=F(n2,1)+f(2,1);
    F(n3,1)=F(n3,1)+f(3,1);
    F(n1+N_n,1)=F(n1+N_n,1)+f(4,1);
    F(n2+N_n,1)=F(n2+N_n,1)+f(5,1);
    F(n3+N_n,1)=F(n3+N_n,1)+f(6,1);
end
%-----------Ensamblaje Matriz Rigidez Elementos de Borde------
syms x
Iso=zeros(6,6); 
Iso(1,1)=double(int((1-x).^2,x,0,1));
Iso(1,2)=double(int((1-x)*x,x,0,1));
Iso(1,3)=0;
Iso(2,1)=double(Iso(1,2));
Iso(2,2)=double(int(x.^2,x,0,1));
Iso(2,3)=0;
Iso(3,1:3)=0;
Iso(1:3,4:6)=Iso(1:3,1:3);
Iso(4:6,1:3)=Iso(1:3,1:3);
Iso(4:6,4:6)=Iso(1:3,1:3);
j=1;
B1=zeros(6,6);
tic
for i=1:1:N_bound
    131
    if Matriz(i,2)==Inf || Matriz(i,2)==-Inf %--Cond. Elementos verticales-
        Matriz(i,2)=0;
        f=2;
    else
    f=1;
    end
    while Matriz(i,1)==elements(j,4)
        n1=elements(j,6);
        n2=elements(j,7);
%------Ensamblaje Matriz Rigidez por Elementos de Borde en cada Elemento--
        B1(1,1)=2*co*((nodes(n2,f)-nodes(n1,f))*Matriz(i,3)*Matriz(i,3)*(1+Matriz(i,2).^2).^0.5)*Iso(1,1);
        B1(1,2)=2*co*((nodes(n2,f)-nodes(n1,f))*Matriz(i,3)*Matriz(i,3)*(1+Matriz(i,2).^2).^0.5)*Iso(1,2);
        B1(1,3)=Iso(1,3);
        B1(1,4)=2*co*((nodes(n2,f)-nodes(n1,f))*Matriz(i,3)*Matriz(i,4)*(1+Matriz(i,2).^2).^0.5)*Iso(1,4);
        B1(1,5)=2*co*((nodes(n2,f)-nodes(n1,f))*Matriz(i,3)*Matriz(i,4)*(1+Matriz(i,2).^2).^0.5)*Iso(1,5);
        B1(1,6)=0;
        B1(2,1)=2*co*((nodes(n2,f)-nodes(n1,f))*Matriz(i,3)*Matriz(i,3)*(1+Matriz(i,2).^2).^0.5)*Iso(2,1);
        B1(2,2)=2*co*((nodes(n2,f)-nodes(n1,f))*Matriz(i,3)*Matriz(i,3)*(1+Matriz(i,2).^2).^0.5)*Iso(2,2);
        B1(2,3)=Iso(1,3);
        B1(2,4)=2*co*((nodes(n2,f)-nodes(n1,f))*Matriz(i,3)*Matriz(i,4)*(1+Matriz(i,2).^2).^0.5)*Iso(2,4);
        B1(2,5)=2*co*((nodes(n2,f)-nodes(n1,f))*Matriz(i,3)*Matriz(i,4)*(1+Matriz(i,2).^2).^0.5)*Iso(2,5);
        B1(2,6)=0;
        B1(3,1:6)=0;
        B1(4,1)=2*co*((nodes(n2,f)-nodes(n1,f))*Matriz(i,3)*Matriz(i,4)*(1+Matriz(i,2).^2).^0.5)*Iso(4,1);
        B1(4,2)=2*co*((nodes(n2,f)-nodes(n1,f))*Matriz(i,3)*Matriz(i,4)*(1+Matriz(i,2).^2).^0.5)*Iso(4,2);
        B1(4,3)=Iso(1,3);
        B1(4,4)=2*co*((nodes(n2,f)-nodes(n1,f))*Matriz(i,4)*Matriz(i,4)*(1+Matriz(i,2).^2).^0.5)*Iso(4,4);
        B1(4,5)=2*co*((nodes(n2,f)-nodes(n1,f))*Matriz(i,4)*Matriz(i,4)*(1+Matriz(i,2).^2).^0.5)*Iso(4,5);
        B1(4,6)=0;
        B1(5,1)=2*co*((nodes(n2,f)-nodes(n1,f))*Matriz(i,3)*Matriz(i,4)*(1+Matriz(i,2).^2).^0.5)*Iso(5,1);
        B1(5,2)=2*co*((nodes(n2,f)-nodes(n1,f))*Matriz(i,3)*Matriz(i,4)*(1+Matriz(i,2).^2).^0.5)*Iso(5,2);
        B1(5,3)=Iso(1,3);
        B1(5,4)=2*co*((nodes(n2,f)-nodes(n1,f))*Matriz(i,4)*Matriz(i,4)*(1+Matriz(i,2).^2).^0.5)*Iso(5,4);
        B1(5,5)=2*co*((nodes(n2,f)-nodes(n1,f))*Matriz(i,4)*Matriz(i,4)*(1+Matriz(i,2).^2).^0.5)*Iso(5,5);
        B1(5,6)=0;
        B1(6,1:6)=0;
        %if n1==6 && n2==84
        % B1
        %end
        %-----------Ensamblaje Matriz de Rigidez Global----------------------------
        K(n1,n1)=K(n1,n1)+B1(1,1);
        K(n1,n2)=K(n1,n2)+B1(1,2);
        K(n1,n1+N_n)=K(n1,n1+N_n)+B1(1,4);
        K(n1,n2+N_n)=K(n1,n2+N_n)+B1(1,5);
        K(n2,n1)=K(n2,n1)+B1(2,1);
        K(n2,n2)=K(n2,n2)+B1(2,2);
        K(n2,n1+N_n)=K(n2,n1+N_n)+B1(2,4);
        K(n2,n2+N_n)=K(n2,n2+N_n)+B1(2,5);
        K(n1+N_n,n1)=K(n1+N_n,n1)+B1(4,1);
        K(n1+N_n,n2)=K(n1+N_n,n2)+B1(4,2);
        K(n1+N_n,n1+N_n)=K(n1+N_n,n1+N_n)+B1(4,4);
        K(n1+N_n,n2+N_n)=K(n1+N_n,n2+N_n)+B1(4,5);
        K(n2+N_n,n1)=K(n2+N_n,n1)+B1(5,1);
        K(n2+N_n,n2)=K(n2+N_n,n2)+B1(5,2);
        K(n2+N_n,n1+N_n)=K(n2+N_n,n1+N_n)+B1(5,4);
        K(n2+N_n,n2+N_n)=K(n2+N_n,n2+N_n)+B1(5,5);
        j=j+1;
    end
end
clear A B1
%--------------------------------------------------------------------------
%---------------VECTOR POTENCIAL MAGNÉTICO ----------------------------
%--------------------------------------------------------------------------
U=sparse(K)\sparse(F);
Ux=U(1:N_n);
Uy=U(N_n+1:2*N_n);
%--------------------------------------------------------------------------
%-------------------GRAFICA FUNCIÓN VECTOR POTENCIAL-----------------------
%--------------------------------------------------------------------------
figure(2)
quiver(two_d_nodes(:,1),two_d_nodes(:,2),Ux,Uy);
title('Vector potencial magnético A');
x=[nodes(3,2) nodes(4,2) nodes(9,2)];
y=[nodes(3,2) nodes(4,2) nodes(9,2)];
%plot(x,y);
U=zeros(N_n,2);
U(1:N_n,1)=Ux;
U(1:N_n,2)=Uy;
clear X Y F1 F2 TF GF g G1 G2 f
%--------------------------------------------------------------------------
%---------------------------PARÁMETRO DE ORDEN-----------------------------
%--------------------------------------------------------------------------
%------Ensamblaje matrices de rigidez por elemento normalizado-------------
syms alf bet 
M1=zeros(3,3);
N(1,1)=1-alf-bet;
N(1,2)=bet;
N(1,3)=alf;
for i=1:1:3
    for j=1:1:3
        M1(i,j)=int(int(N(1,i)*N(1,j),alf,0,1-bet),bet,0,1);
    end
end
G11=Ks*Ks*double(int(int(transpose(N(1,:))*N(1,:)*N(1,1)*N(1,1),alf,0,1-bet),bet,0,1));
G12=Ks*Ks*2*double(int(int(transpose(N(1,:))*N(1,:)*N(1,1)*N(1,2),alf,0,1-bet),bet,0,1));
G13=Ks*Ks*2*double(int(int(transpose(N(1,:))*N(1,:)*N(1,1)*N(1,3),alf,0,1-bet),bet,0,1));
G14=Ks*Ks*double(int(int(transpose(N(1,:))*N(1,:)*N(1,3)*N(1,3),alf,0,1-bet),bet,0,1));
G15=Ks*Ks*double(int(int(transpose(N(1,:))*N(1,:)*N(1,2)*N(1,2),alf,0,1-bet),bet,0,1));
G16=Ks*Ks*2*double(int(int(transpose(N(1,:))*N(1,:)*N(1,2)*N(1,3),alf,0,1-bet),bet,0,1));
F11=zeros(1,3);
F12=zeros(1,3);
F13=zeros(1,3);
F14=zeros(3,1);
F15=zeros(3,1);
F16=zeros(3,1);
for i=1:1:1
    for j=1:1:3
        F11(i,j)=Ks*int(int(N(1,1)*N(1,j),alf,0,1-bet),bet,0,1);
        F12(i,j)=Ks*int(int(N(1,2)*N(1,j),alf,0,1-bet),bet,0,1);
        F13(i,j)=Ks*int(int(N(1,3)*N(1,j),alf,0,1-bet),bet,0,1);
        F14(j,1)=Ks*int(int(N(1,1)*N(1,j),alf,0,1-bet),bet,0,1);
        F15(j,1)=Ks*int(int(N(1,2)*N(1,j),alf,0,1-bet),bet,0,1);
        F16(j,1)=Ks*int(int(N(1,3)*N(1,j),alf,0,1-bet),bet,0,1);
    end
end
%------------Ensamblaje Matriz de rigidez global---------------------------
clear C0 C1 C C2 A d
d = ones(3,3);
C0=zeros(3,2);
C1=zeros(3,3);
C2=zeros(3,3);
C=zeros(3,3);
A=zeros(6,6);
Kd=sparse(N_n*2,N_n*2);
Ki=sparse(N_n*2,N_n*2);
for i = 1:1:N_e-one_d_elements
    n1 = two_d_elements(i,1); n2 = two_d_elements(i,2); n3 = two_d_elements(i,3);
    d(1,2) = nodes(n1,1); d(1,3) = nodes(n1,2);
    d(2,2) = nodes(n2,1); d(2,3) = nodes(n2,2);
    d(3,2) = nodes(n3,1); d(3,3) = nodes(n3,2);
    J =(det(d));
    %-----Jacobiano por elemento-------------------
    %-----------Ensamblaje matriz A-------------------------------------------
    C0(1,1)=(nodes(n2,2)-nodes(n3,2))/J;
    C0(1,2)=(nodes(n3,1)-nodes(n2,1))/J;
    C0(2,1)=(nodes(n3,2)-nodes(n1,2))/J;
    C0(2,2)=(nodes(n1,1)-nodes(n3,1))/J;
    C0(3,1)=(nodes(n1,2)-nodes(n2,2))/J;
    C0(3,2)=(nodes(n2,1)-nodes(n1,1))/J;
    C1=C0*transpose(C0)*J*0.5;
    C2(1:3,1:3)=G11*(U(n1,1).^2+U(n1,2).^2)+G12*(U(n1,1)*U(n2,1)+U(n1,2)*U(n2,2)...
    )+G13*(U(n1,1)*U(n3,1)+U(n1,2)*U(n3,2))+G14*(U(n3,1).^2+U(n3,2).^2)+G16*(U...
    (n2,1)*U(n3,1)+U(n2,2)*U(n3,2))+G15*(U(n2,1).^2+U(n2,2).^2);
    C=C1+C2*J;
    A(1:3,1:3)=C;
    A(4:6,4:6)=C;
    A(1:3,4:6)=U(n1,1)*C0(:,1)*F11*J+U(n1,2)*C0(:,2)*F11*J+U(n2,1)*C0(:,1)*F12*J+...
    U(n2,2)*C0(:,2)*F12*J+U(n3,1)*C0(:,1)*F13*J+U(n3,2)*C0(:,2)*F13*J-...
    U(n1,1)*F14*transpose(C0(:,1))*J-U(n1,2)*F14*transpose(C0(:,2))*J-...
    U(n2,1)*F15*transpose(C0(:,1))*J-U(n2,2)*F15*transpose(C0(:,2))*J-...
    U(n3,1)*F16*transpose(C0(:,1))*J-U(n3,2)*F16*transpose(C0(:,2))*J;...
    A(4:6,1:3)=transpose(A(1:3,4:6));
    %------------Ensamblaje matriz Ki y Kd------------------------------------
    nump=N_n;
    Ki(n1,n1)=Ki(n1,n1)+A(1,1);
    Ki(n1,n2)=Ki(n1,n2)+A(1,2);
    Ki(n1,n3)=Ki(n1,n3)+A(1,3);
    Ki(n2,n1)=Ki(n2,n1)+A(2,1);
    Ki(n2,n2)=Ki(n2,n2)+A(2,2);
    Ki(n2,n3)=Ki(n2,n3)+A(2,3);
    Ki(n3,n1)=Ki(n3,n1)+A(3,1);
    Ki(n3,n2)=Ki(n3,n2)+A(3,2);
    Ki(n3,n3)=Ki(n3,n3)+A(3,3);
    Ki(n1,n1+N_n)=Ki(n1,n1+N_n)+A(1,4);
    Ki(n1,n2+N_n)=Ki(n1,n2+N_n)+A(1,5);
    Ki(n1,n3+N_n)=Ki(n1,n3+N_n)+A(1,6);
    Ki(n2,n1+N_n)=Ki(n2,n1+N_n)+A(2,4);
    Ki(n2,n2+N_n)=Ki(n2,n2+N_n)+A(2,5);
    Ki(n2,n3+N_n)=Ki(n2,n3+N_n)+A(2,6);
    Ki(n3,n1+N_n)=Ki(n3,n1+N_n)+A(3,4);
    Ki(n3,n2+N_n)=Ki(n3,n2+N_n)+A(3,5);
    Ki(n3,n3+N_n)=Ki(n3,n3+N_n)+A(3,6);
    Ki(n1+nump,n1)= Ki(n1+nump,n1)+A(4,1);
    Ki(n1+nump,n2)= Ki(n1+nump,n2)+A(4,2);
    Ki(n1+nump,n3)= Ki(n1+nump,n3)+A(4,3);
    Ki(n1+nump,n1+nump)= Ki(n1+nump,n1+nump)+A(4,4);
    Ki(n1+nump,n2+nump)= Ki(n1+nump,n2+nump)+A(4,5);
    Ki(n1+nump,n3+nump)= Ki(n1+nump,n3+nump)+A(4,6);
    Ki(n2+nump,n1)= Ki(n2+nump,n1)+A(5,1);
    Ki(n2+nump,n2)= Ki(n2+nump,n2)+A(5,2);
    Ki(n2+nump,n3)= Ki(n2+nump,n3)+A(5,3);
    Ki(n2+nump,n1+nump)= Ki(n2+nump,n1+nump)+A(5,4);
    Ki(n2+nump,n2+nump)= Ki(n2+nump,n2+nump)+A(5,5);
    Ki(n2+nump,n3+nump)= Ki(n2+nump,n3+nump)+A(5,6);
    Ki(n3+nump,n1)= Ki(n3+nump,n1)+A(6,1);
    Ki(n3+nump,n2)= Ki(n3+nump,n2)+A(6,2);
    Ki(n3+nump,n3)= Ki(n3+nump,n3)+A(6,3);
    Ki(n3+nump,n1+nump)= Ki(n3+nump,n1+nump)+A(6,4);
    Ki(n3+nump,n2+nump)= Ki(n3+nump,n2+nump)+A(6,5);
    Ki(n3+nump,n3+nump)= Ki(n3+nump,n3+nump)+A(6,6);
    Kd(n1,n1)=Kd(n1,n1)+J*M1(1,1);
    Kd(n1,n2)=Kd(n1,n2)+J*M1(1,2);
    Kd(n1,n3)=Kd(n1,n3)+J*M1(1,3);
    Kd(n2,n1)=Kd(n2,n1)+J*M1(2,1);
    Kd(n2,n2)=Kd(n2,n2)+J*M1(2,2);
    Kd(n2,n3)=Kd(n2,n3)+J*M1(2,3);
    Kd(n3,n1)=Kd(n3,n1)+J*M1(3,1);
    Kd(n3,n2)=Kd(n3,n2)+J*M1(3,2);
    Kd(n3,n3)=Kd(n3,n3)+J*M1(3,3);
    Kd(n1,n1+nump)=Kd(n1,n1+nump);
    Kd(n1,n2+nump)=Kd(n1,n2+nump);
    Kd(n1,n3+nump)=Kd(n1,n3+nump);
    Kd(n2,n1+nump)=Kd(n2,n1+nump);
    Kd(n2,n2+nump)=Kd(n2,n2+nump);
    Kd(n2,n3+nump)=Kd(n2,n3+nump);
    Kd(n3,n1+nump)=Kd(n3,n1+nump);
    Kd(n3,n2+nump)=Kd(n3,n2+nump);
    Kd(n3,n3+nump)=Kd(n3,n3+nump);
    Kd(n1+nump,n1)= Kd(n1+nump,n1);
    Kd(n1+nump,n2)= Kd(n1+nump,n2);
    Kd(n1+nump,n3)= Kd(n1+nump,n3);
    Kd(n1+nump,n1+nump)= Kd(n1+nump,n1+nump)+J*M1(1,1);
    Kd(n1+nump,n2+nump)= Kd(n1+nump,n2+nump)+J*M1(1,2);
    Kd(n1+nump,n3+nump)= Kd(n1+nump,n3+nump)+J*M1(1,3);
    Kd(n2+nump,n1)= Kd(n2+nump,n1);
    Kd(n2+nump,n2)= Kd(n2+nump,n2);
    Kd(n2+nump,n3)= Kd(n2+nump,n3);
    Kd(n2+nump,n1+nump)= Kd(n2+nump,n1+nump)+J*M1(2,1);
    Kd(n2+nump,n2+nump)= Kd(n2+nump,n2+nump)+J*M1(2,2);
    Kd(n2+nump,n3+nump)= Kd(n2+nump,n3+nump)+J*M1(2,3);
    Kd(n3+nump,n1)= Kd(n3+nump,n1);
    Kd(n3+nump,n2)= Kd(n3+nump,n2);
    Kd(n3+nump,n3)= Kd(n3+nump,n3);
    Kd(n3+nump,n1+nump)= Kd(n3+nump,n1+nump)+J*M1(3,1);
    Kd(n3+nump,n2+nump)= Kd(n3+nump,n2+nump)+J*M1(3,2);
    Kd(n3+nump,n3+nump)= Kd(n3+nump,n3+nump)+J*M1(3,3);
    clear A
end
%--------------------------------------------------------------------------
%-------------------CÁLCULO DE EIGENVALORES Y EIGENVECTORES------------
%--------------------------------------------------------------------------
[V2,D2]=eigs(Ki,Kd,10,0);
D3=zeros(1,10);
for i=1:1:10
    D3(1,i)=D2(i,i);
end
a11=min(D3);
for i=1:1:10
    if D2(i,i)==a11
        a12=i;
    end
end
a11=a12;
toc
%--------------------------------------------------------------------------
%------------------GRÁFICA DEL PARÁMETRO DE ORDEN--------------------------
%--------------------------------------------------------------------------
Mo1=sparse(nump,2*nump);
for i=1:1:nump
    Mo1(i,i)=V2(i,a11);
    Mo1(i,i+nump)=V2(i+nump,a11);
end
Mo=Mo1*V2(:,a11);
%--------------Grafica parámetro de orden. Superficie----------------------
figure(3)
trisurf(two_d_elements,two_d_nodes(:,1),two_d_nodes(:,2),Mo.^0.5,'Edgecolor','none')
shading interp
title('Parametro de Orden \Psi GL','fontsize',14)
%--------------------------------------------------------------------------
%---------------------FINAL CÓDIGO-----------------------------------
