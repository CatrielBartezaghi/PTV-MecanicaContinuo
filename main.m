clc;
clear all;
%% Datos del problema
grPol = 3;
v=100;  
%
modVelBusc = 33;
posVelBusc = [16,12,5];

%% Datos gráfico
%Número de divisiones.
nPartx = 10;
nParty = 10;
nPartz = 10;
%Número de divisiones de la línea O-A.
nPartOA = (nPartx+nParty+nPartz)/2-4;
%% Matriz simbólica de producto coordenadas polinomio
syms x y z

nCoef=(grPol+1)*(grPol+2)*(grPol+3)/6;
ms_N = sym(zeros(nCoef,1));
ind = 0;
for iz = 0:grPol
 for iy = 0:grPol-iz
   for ix = 0:grPol-iz-iy
    ind = ind+1;
    ms_N(ind) = (x^2+y^2-400)*x^ix*y^iy*z^iz;
   end
 end
end
%% Matriz de deformación simbólica
syms ms_B
ms_B = simplify([diff(ms_N.',x);diff(ms_N.',y);diff(ms_N.',z)]);
%% Matriz de rigidez (matriz de coeficiente del sistema lineal)
syms ms_K

funcion_arriba=sqrt(400-(x^2));
ms_K = int(int(int(ms_B.'*ms_B,y,0,funcion_arriba),x,0,20),z,0,10);
m_K = double(ms_K);
%vpa(ms_K,10);

%% Vector de fuerza (vector de términos independientes)
syms ms_Fr ms_F

ms_Fr = subs(ms_N,[x,y,z],[0,0,10]);
ms_F = sym(v)*ms_Fr;
m_F = double(ms_F);
%m_F = vpa(ms_F);
%% Resolución del sistema en forma simbólica
syms ms_A
ms_A = ms_K\ms_F;
%% Resolución del sistema en forma numérica
m_A = m_K\m_F;
%% Gráfico de la función potencial
figure(1)

%La grilla tiene que ser 3d.
[m_CoordX,m_CoordY,m_CoordZ] = meshgrid(linspace(0,20,2*nPartx+1),linspace(0,20,2*nParty+1),linspace(0,10,2*nPartz+1));
%Para graficar la L se pone NaN en las coordenadas fuera del dominio.
m_CoordYAux=m_CoordY;
m_CoordXAux=m_CoordX;
m_CoordZAux=m_CoordZ;
m_CoordY(m_CoordYAux>sqrt(400-(m_CoordXAux.^2))) = NaN;
m_CoordY(m_CoordYAux<-sqrt(400-(m_CoordXAux.^2))) = NaN;
m_CoordY(m_CoordYAux<0) = NaN;
m_CoordX(m_CoordXAux<0) = NaN;

nPuntos = numel(m_CoordX);

m_XGraf = zeros(nPuntos,nCoef);
%Se evalua cada uno de los elementos de ms_X para todos los puntos de grilla de la gráfica.
for iCoef = 1:nCoef
 m_XGraf(:,iCoef) = subs(ms_N(iCoef),{x,y,z},{m_CoordX(:),m_CoordY(:),m_CoordZ(:)});
 
end
%Función potencial
m_FunPot = reshape(m_XGraf*(m_A),size(m_CoordX));
%Gráfico
for i=1:size(m_CoordX(:,1,1))
    figure(i)
    mesh(m_CoordX(:,:,i),m_CoordY(:,:,i),m_CoordZ(:,:,i),m_FunPot(:,:,i)) 
    xlabel('Eje X'); ylabel('Eje Y'); zlabel('Eje Z');
end


%% Gráfico de módulo de velocidad sobre la línea OA
hold on
figure()
hold on
m_CoordXOA = linspace(0,0,nPartOA);
m_CoordYOA = linspace(0,0,nPartOA);
m_CoordZOA = linspace(0,10,nPartOA);
nPuntos = numel(m_CoordXOA);
%Se divide en dos líneas la matriz de deformación.
m_BGraf1 = zeros(nPuntos,nCoef);
m_BGraf2 = zeros(nPuntos,nCoef);
m_BGraf3 = zeros(nPuntos,nCoef);
%Se evalua cada uno de los elementos de ms_B para todos los puntos degrilla de la gráfica.
for iCoef = 1:nCoef
 m_BGraf1(:,iCoef) = subs(ms_B(1,iCoef),{x,y,z},{m_CoordXOA,m_CoordYOA,m_CoordZOA});
 m_BGraf2(:,iCoef) = subs(ms_B(2,iCoef),{x,y,z},{m_CoordXOA,m_CoordYOA,m_CoordZOA});
 m_BGraf3(:,iCoef) = subs(ms_B(3,iCoef),{x,y,z},{m_CoordXOA,m_CoordYOA,m_CoordZOA});
end
%Vector de velocidad
m_VecVelX = -m_BGraf1*m_A;
m_VecVelY = -m_BGraf2*m_A;
m_VecVelZ = -m_BGraf3*m_A;
%Gráfico del módulo de velocidad
%Se grafica en función de las coordenadas x de los puntos.

mod = sqrt(abs(m_VecVelX).^2+abs(m_VecVelY).^2+abs(m_VecVelZ).^2);
plot(m_CoordZOA,mod,'LineWidth',1.5);

xlabel('Coordenada Z'); ylabel('Módulo de la velocidad');
 
 
 %% Determinación de la velocidad de entrada para que un cierto punto tenga un módulo de velocidad dado
modVelBusc = 33;
posVelBusc = [16,12,5];
m_BGrafposVelBusc = subs(ms_B,{x,y,z},{posVelBusc(1),posVelBusc(2),posVelBusc(3)});
m_VecVelposVelBusc = m_BGrafposVelBusc*m_A;
absVelEntr = modVelBusc/norm(m_VecVelposVelBusc)*abs(v)


%% Gráfico de la distribución de vectores de velocidad en el dominio

figure();
%Datos para graficar el cuarto de circunferencia
t = 0:pi/30:pi/2;
xc = 20*cos(t);
yc = 20*sin(t);
z0 = zeros(1,16);
z10 = 10*ones(1,16);
%Gráfico del contorno de la figura
plot3([0,20,20 ,0,0 ,0,0,0],...
      [0,0 ,0 ,0 ,20,20,0,0],...
      [0,0 ,10,10,10 ,0,0,10],'k','LineWidth',2)
hold on
plot3(xc,yc,z0,'k','LineWidth',2);
plot3(xc,yc,z10,'k','LineWidth',2);
axis square

[m_CoordX,m_CoordY,m_CoordZ] = meshgrid(linspace(0,20,11),linspace(0,20,11),linspace(0,10,11));
m_CoordYAux=m_CoordY;
m_CoordXAux=m_CoordX;
m_CoordZAux=m_CoordZ;
m_CoordY(m_CoordYAux>sqrt(400-(m_CoordXAux.^2))) = NaN;
m_CoordY(m_CoordYAux<-sqrt(400-(m_CoordXAux.^2))) = NaN;
m_CoordY(m_CoordYAux<0) = NaN;
m_CoordX(m_CoordXAux<0) = NaN;

nPuntos = numel(m_CoordX);
%Se divide en dos líneas la matriz de gradientes
m_BGraf1 = zeros(nPuntos,nCoef);
m_BGraf2 = zeros(nPuntos,nCoef);
m_BGraf3 = zeros(nPuntos,nCoef);
%Se evalua cada uno de los elementos de ms_B para todos los puntos de grilla de la gráfica.
for iCoef = 1:nCoef
m_BGraf1(:,iCoef) = subs(ms_B(1,iCoef),{x,y,z},{m_CoordX(:),m_CoordY(:),m_CoordZ(:)});
m_BGraf2(:,iCoef) = subs(ms_B(2,iCoef),{x,y,z},{m_CoordX(:),m_CoordY(:),m_CoordZ(:)});
m_BGraf3(:,iCoef) = subs(ms_B(3,iCoef),{x,y,z},{m_CoordX(:),m_CoordY(:),m_CoordZ(:)});
end
%Vector de flujo
m_VecVelX = -m_BGraf1*m_A;
m_VecVelY = -m_BGraf2*m_A;
m_VecVelZ = -m_BGraf3*m_A;

%El quiver escala a partir del tamaño automático que pone automáticamente
quiver3(m_CoordX(:),m_CoordY(:),m_CoordZ(:),m_VecVelX,m_VecVelY,m_VecVelZ,5);
hold off
axis equal


