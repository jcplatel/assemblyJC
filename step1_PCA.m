function [Mnew,PCout] = step1_PCA(M,Npc_max)
%STEP1_PCA Apply Principal Component Analysis
%   INPUTS
%   M               matrix Npix X Nt - movie image
%   Npc_max     	scalar - number max of principal components
%   OUTPUTS
%   Mnew            matrix Npix X Nt - centralized movie image used for PCA
%   PCout           struct with fields:
%                   - pc_XFilter (Npix X Npc)   spatial eigenvectors
%                   - pc_TimeCourse (Nt X Npc)  temporal eigenvectors
%                   - pc_EigVal (Npc)           eigenvalues


% -- Number of pixel and time frames
[Npix Nt]= size(M);


% -- Centralize data
fluo_pix_mt = mean(M,2); % average fluo intensity for each pixel over time
M_mt = repmat(fluo_pix_mt,1,Nt);
Mn_1 = (M - M_mt); % centralize (temporal)
clear M M_mt

fluo_t_mpix = mean(Mn_1,1); % average fluo intensity at each time frame - average over all pixels
M_mpix = repmat(fluo_t_mpix,Npix,1);
Mn_2 = Mn_1 - M_mpix; % centralize (spatial)
clear Mn_1 M_mpix

Mnew = Mn_2; % final data Npix X Nt
clear Mn_2

 
% -- Singular Value Decomposition  
% M=UDV'  eigenvectors U: spatial filters    V: signal time courses
% Calculate eigenvectors on smaller dimension, then compute the
% eigenvectors in the complementary space (temporal, spatial)
% U and V are both orthonormal

if Nt<=Npix             % spatial covariance on temporal matrix
    
    % find V
%     MtM = (Mnew'*Mnew)/Npix;
    MtM = (Mnew(:,1:end-1)'*Mnew(:,2:end)+Mnew(:,2:end)'*Mnew(:,1:end-1))/Npix/2;
    % eigenvectors and eigenvalues
    [V_all,EigVal_all] = eigs(MtM,Npc_max);
    % keep positive non-zero eig val
    [i_pos,j_pos] = find(EigVal_all > eps);
    EigVal = EigVal_all(i_pos,j_pos);
    % keep corresponding eigenvectors
    V = V_all(:,j_pos);
   
    % find U
    mat_SVD = (Npix*EigVal).^(1/2);
    U =  Mnew(:,1:end-1) * V / mat_SVD;
    
else                    % temporal covariance on spatial matrix
    
    % find U
%     MMt = (Mnew*Mnew')/Nt;
    MMt = (Mnew(:,1:end-1)*Mnew(:,2:end)'+Mnew(:,2:end)*Mnew(:,1:end-1)')/Nt/2;
    % eigenvectors and eigenvalues
    [U_all,EigVal_all] = eigs(MMt,Npc_max);
    % keep positive non-zero eig val
    [i_pos,j_pos] = find(EigVal_all > eps);
    EigVal = EigVal_all(i_pos,j_pos);
    % keep corresponding eigenvectors
    U = U_all(:,j_pos);
   
    % find V
    mat_SVD = (Nt*EigVal).^(1/2);
    V = (mat_SVD\U'*Mnew)';
    
end



% -- OUTPUT
PCout.pc_XFilter = U;
PCout.pc_TimeCourse = V;
PCout.pc_EigVal = diag(EigVal);


end



