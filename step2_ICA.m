function [ ICout ] = step2_ICA(STmu,UnmixingMat_init,pc_XFilter,pc_TimeCourse,Nsrc,Niter_max,TermTol)
%STEP2_ICA Apply Independent Component Analysis
%   INPUTS
%   STmu            scalar - relative coefficient for spatio-temporal ICA
%                   (0=temporal, 1=spatial)
%   UnmixingMat_init 
%                   matrix Nic X Npc - initialization of unmixing matrix
%   pc_XFilter      matrix Npix X Npc - spatial eigenvectors
%   pc_TimeCourse   matrix Nt X Npc - temporal eigenvectors
%   Nsrc            scalar - number max of independent components
%   Niter_max       scalar - number max of iteration steps
%   TermTol         scalar - termination tolerance
%   OUTPUTS
%   ICout           struct with fields:
%                   - ic_XFilter (Npix X Nic)   spatial ic
%                   - skw_spatial (Nic)         spatial skewness
%                   - ic_TimeCourse (Nt X Nic)  temporal ic
%                   - skw_temporal (Nic)        temporal skewness
%                   - UnmixingMat (Nic X Npc)   unmixing matrix solution
%                   - Niter                     number of steps for cv
%                   - UnmixingMat_err           unmixing matrix final error


% Note: Xfilter and TimeCourse pc must have been selected


% -- Number of principal components
[~,Npc] = size(pc_TimeCourse);


% -- Check inputs
if isempty(Nsrc)
    Nsrc = Npc;
end
if isempty(Niter_max)
    Niter_max = 100;
end
if isempty(TermTol) || nargin<7
    TermTol = 1e-6;
end
if isempty(STmu)
    STmu = 0;
end
if size(UnmixingMat_init,1)~=Nsrc || size(UnmixingMat_init,2)~=Npc
    error('Initialization unmixing matrix -> wrong size')
end
if Npc<Nsrc
    error('More ICs than PCs')
end


% -- Select type of ICA
if STmu == 0  % independent signal time course sources
    
    mixed_sig = pc_TimeCourse; % mixed signal (Nt X Npc)
    
    fprintf('temporal ICA \n')

elseif STmu == 1        % independent spatial sources
    
    mixed_sig = pc_XFilter; % mixed signal (Npix X Npc)
    
    fprintf('spatial ICA \n')
    
else % spatio-temporal ICA
    
    mixed_sig = [STmu*pc_XFilter ;(1-STmu)*pc_TimeCourse]; % spatio-temporal ( Npix+Nt X Npc)
    
    fprintf('spatio-temporal ICA mu = %d \n',STmu)
    
    
end

% -- Dimension of mixed signal
Ni = size(mixed_sig,1);


% -- Compute ICA
mixed_sig = mixed_sig - repmat(mean(mixed_sig,1),[Ni 1]); % center each pc
[UnmixingMat,Niter,UnmixingMat_err] = compute_ICA();


% -- Determine sources (unmixed signals)
ic_Sources = mixed_sig*UnmixingMat'; % (Ni[Nt or Npix or Nt+Npix] X Nsrc)


% -- Evaluate skewness of each source
skw_unmixed = skewness(ic_Sources,0);


% -- Retrieve independent components
if STmu == 0            % temporal ICA
    
    ic_TimeCourse = ic_Sources;
    ic_XFilter = pc_XFilter*UnmixingMat';
    skw_temporal = skw_unmixed;
    skw_spatial = skewness(ic_XFilter,0);
    
elseif STmu == 1        % spatial ICA
    
    ic_XFilter = ic_Sources;
    ic_TimeCourse = pc_TimeCourse*UnmixingMat';
    skw_temporal = skewness(ic_TimeCourse,0);
    skw_spatial = skw_unmixed;
    
else                    % spatio-temporal ICA
  
    ic_XFilter = ic_Sources(1:size(pc_XFilter,1),:);
    ic_TimeCourse = ic_Sources(size(pc_XFilter,1)+1:end,:);
    skw_temporal = skewness(ic_TimeCourse,0);
    skw_spatial = skewness(ic_XFilter,0);
    
    
end


% -- Sort sources by temporal skewness
[skw_temporal,IX] = sort(skw_temporal, 'descend');
ic_XFilter = ic_XFilter(:,IX);
ic_TimeCourse = ic_TimeCourse(:,IX);
skw_spatial = skw_spatial(:,IX);



% -- Nested called functions

    function [UnmixingMat,Niter,UnmixingMat_err] = compute_ICA()
        
        UnmixingMat = UnmixingMat_init;
        
        % orthogonalize
        UnmixingMat = (UnmixingMat*UnmixingMat')^(1/2)\UnmixingMat;
        
        i_iter = 0;
        UnmixingMat_err = zeros(Niter_max,1);
        cv_test = 100;
        
        while ( (i_iter < Niter_max) && (cv_test > TermTol) )
            
            i_iter = i_iter+1;
            
            % save matrix iteration n-1
            UnmixingMat_old = UnmixingMat;
 
            % estimate source
            unmixed_sig = mixed_sig*UnmixingMat'; % (Ni[Nt or Npix] X Nsrc)
    
            unmixed_sig2 = unmixed_sig'.^2;
            
            % update unmixing matrix
            UnmixingMat = unmixed_sig2 * mixed_sig;
            
            % orthogonalize
            UnmixingMat = (UnmixingMat*UnmixingMat')^(1/2)\UnmixingMat;
            
            % test for termination convergence
            cv_test = 1 - min(abs(diag(UnmixingMat*UnmixingMat_old')));
            
            UnmixingMat_err(i_iter) = cv_test;
            
        end
        
        Niter = i_iter;
        
        if Niter < Niter_max
            fprintf('Convergence in %d iterations\n', Niter)
        else
            fprintf('Failed to converge in %d iterations; current error %3.3g \n', ...
                Niter, cv_test)
        end
        
    end


% -- OUTPUT
ICout.ic_XFilter = ic_XFilter;
ICout.skw_spatial = skw_spatial;
ICout.ic_TimeCourse = ic_TimeCourse;
ICout.skw_temporal = skw_temporal;
ICout.UnmixingMat = UnmixingMat;
ICout.Niter = Niter;
ICout.UnmixingMat_err = UnmixingMat_err;


end