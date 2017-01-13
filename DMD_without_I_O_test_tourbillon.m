clear all
timestep = 250*3.5e-9;
diameter = 3.78*1.67e-3;
speed = 937.7594;
nombre_total_fichiers=256;
u=10;
M_u = u*ones(200,250);
M_v = u*ones(200,250);
Omega = 2;

for i=51:150,
    for j=1:100,
        M_v(i,j) = u+Omega*(j-50)/sqrt(i*i+j*j);
        M_u(i,j) = u+Omega*(i-100)/sqrt(i*i+j*j);
    end
end

M_U = M_u(:,1);
M_V = M_v(:,1);

for i=2:size(M_u,2),
    M_U = vertcat(M_U, M_u(:,i));
    M_V = vertcat(M_V, M_v(:,i));
end

% LECTURE DE L'INTEGRALITE DU REPERTOIRE

for i=1:nombre_total_fichiers,
    
    
    for j=2:size(M_u,2),
        
        if j>i,
            M_v(:,j) = M_v(:,j-1);
            M_u(:,j) = M_u(:,j-1);
        end
    
        if j<i,
            M_v(:,j) = u;
            M_u(:,j) = u;
        end
    
        if j==size(M_u,2),
            M_U_temp = M_u(:,1);
            M_V_temp = M_v(:,1);
            for k=2:size(M_u,2),
                M_U_temp = vertcat(M_U_temp, M_u(:,k));
                M_V_temp = vertcat(M_V_temp, M_v(:,k));
            end
        end
    end
    
    if i < nombre_total_fichiers,
        U1(:,i)=M_U_temp;
        V1(:,i)=M_V_temp;
    end
    
    if i > 1,
        U2(:,i-1)=M_U_temp;
        V2(:,i-1)=M_V_temp;
    end
    
end

% calcul sur les champs fluctuants
% U1mean=mean(U1,2);
% V1mean=mean(V1,2);
% U2mean=mean(U2,2);
% V2mean=mean(V2,2);
% 
% for i=1:nombre_total_fichiers-1,
%     U1(:,i)=U1(:,i)-U1mean;
%     V1(:,i)=V1(:,i)-V1mean;
%     U2(:,i)=U2(:,i)-U2mean;
%     V2(:,i)=V2(:,i)-V2mean;
% end
    

% RESOLUTION DU SYSTEME LINEAIRE
Uop = linsolve(U1,U2);
Vop = linsolve(V1,V2);

% CALCUL DES VALEURS ET VECTEURS PROPRES
[U_vector,U_lambda] = eig(Uop);
[V_vector,V_lambda] = eig(Vop);

% CALCUL DES MODES DMD
% calcul "boucle"
%Udmd(size(X,1),nombre_total_fichiers-1) = 0;
%Vdmd(size(X,1),nombre_total_fichiers-1) = 0;
%for i = 1:size(X,1),
%    for j = 1:nombre_total_fichiers-1,
%         for k = 1:nombre_total_fichiers-1,
%             Udmd(i,j)=Udmd(i,j)+U_vector(k,j).*U1(i,k);
%             Vdmd(i,j)=Vdmd(i,j)+V_vector(k,j).*V1(i,k);
%         end
%    end
%end
% calcul matriciel
Udmd = U1*U_vector;
Vdmd = V1*V_vector;

% CALCUL DE L'INTENSITE DES MODES
for i = 1:nombre_total_fichiers-1,
    U_mode_energy(i) = norm(Udmd(:,i));
    V_mode_energy(i) = norm(Vdmd(:,i));
end
 
% CLASSEMENT DES MODES (suivant leur intensité)
[U_mode_energy, Uindices] = sort(U_mode_energy, 'descend');
[V_mode_energy, Vindices] = sort(V_mode_energy, 'descend');
Utemp = Udmd(:,Uindices);
Vtemp = Vdmd(:,Vindices);
Udmd = Utemp;
Vdmd = Vtemp;
clear Utemp Vtemp;
Utemp = U_vector(:,Uindices);
Vtemp = V_vector(:,Vindices);
U_vector = Utemp;
V_vector = Vtemp;
clear Utemp Vtemp;

% CLASSEMENT DES VALEURS PROPRES
for i=1:nombre_total_fichiers-1,
    Utemp(i,i) = U_lambda(Uindices(i),Uindices(i));
    Vtemp(i,i) = V_lambda(Vindices(i),Vindices(i));
end
U_lambda = Utemp;
V_lambda = Vtemp;
clear Utemp Vtemp;  

% CALCUL DES FREQUENCES
for i = 1:nombre_total_fichiers-1,
    U_mode_frequency(i) = imag(log(U_lambda(i,i)))/timestep/2/pi*diameter/speed;
    V_mode_frequency(i) = imag(log(V_lambda(i,i)))/timestep/2/pi*diameter/speed;
end

% CALCUL DES COEFFICIENTS TEMPORELS DES DIFFERENTS MODES
for i=1:nombre_total_fichiers-1,
    Uphi(:,i)=Udmd(:,i)./norm(Udmd(:,i));
    Vphi(:,i)=Vdmd(:,i)./norm(Vdmd(:,i));
end
Ua=U1'*Uphi;
Va=V1'*Vphi;

% SEPARATION PARTIES REELLES ET COMPLEXES
Ureal = real(Udmd);
Uimag = imag(Udmd);
Vreal = real(Vdmd);
Vimag = imag(Vdmd);

% % ecriture des valeurs propres
% fid = fopen(['LAMBDA.dat'],'w');
% for i = 1:nombre_total_fichiers-1,
%     fprintf(fid,'%f %f %f %f %f %f %f %f %f %f\n',real(U_lambda(i,i)), imag(U_lambda(i,i)),...
%         real(log(U_lambda(i,i)))/timestep, imag(log(U_lambda(i,i)))/timestep, U_mode_energy(i), ...
%         real(V_lambda(i,i)), imag(V_lambda(i,i)),...
%         real(log(V_lambda(i,i)))/timestep, imag(log(V_lambda(i,i)))/timestep, V_mode_energy(i));
% end
% fclose(fid);
% % ecriture des fréquences et amplitudes
% fid = fopen(['INTENSITE.dat'],'w');
% for i = 1:nombre_total_fichiers-1,
%     fprintf(fid,'%f %f %f %f\n',U_mode_frequency(i), U_mode_energy(i),V_mode_frequency(i), V_mode_energy(i));
% end
% fclose(fid);

