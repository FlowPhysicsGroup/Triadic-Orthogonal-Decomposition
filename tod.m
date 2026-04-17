function [L,P,f,T,A,Xi] = tod(X,varargin)
% TOD triadic orthogonal decomposition
%
%            fn
%             ^
%    fk_      |________
%     |\     /|       |
%        \ /  |       |
%        / \  |       |
%      /     \|       |
% ----+-------+-------+-> fl
%     |       |\     / 
%     |       |  \ /   
%     |       |  / \  
%     |_______|/     \
%             |
% Figure. Sketch of the fl - fn plane. Triads are expressed as frequency
%         triplets, {fk=fn-l,fl,fn}, or frequency index triplets,
%         (k=n-l,l,n). The principal region consists of fn>=0 for real
%         data, all triads otherwise.
%
%  [L,P,F] = TOD(X) returns the triadic orthogonal decomposition of the
%  data matrix X. The first dimension of X must be time, the second
%  dimension variable indices. X can have any number of additional spatial
%  dimensions. The mode bispectrum is returned in L and the modes in cell
%  array P. The spatial dimensions of the modes are identical to those of
%  X. The convective modes are stored in P{l,n}(1,...) and the recipient
%  modes in P{l,n}(2,...). The cell indices {l,n} of P are the donor and
%  recipient indices. F is the frequency vector. If DT is not specified,
%  the frequency index is returned in F. Although TOD(X) automatically
%  chooses default spectral estimation parameters, it is recommended to
%  manually specify problem-dependent parameters on a case-to-case basis.
%
%  [L,P,F] = TOD(X,WINDOW) uses a temporal window. If WINDOW is a
%  vector, X is divided into segments of the same length as WINDOW. Each
%  segment is then weighted (pointwise multiplied) by WINDOW. If WINDOW is
%  a scalar, a Hamming window of length WINDOW is used. If WINDOW is
%  omitted or empty, a Hamming window is used.
%
%  [L,P,F] = TOD(X,WINDOW,WEIGHT) uses a spatial inner product weight,
%  usually quadrature weights. WEIGHT must have the same spatial dimensions
%  as X.
%
%  [L,P,F] = TOD(X,WINDOW,WEIGHT,NOVERLAP) increases the number of
%  segments by overlapping consecutive blocks by NOVERLAP snapshots.
%  NOVERLAP defaults to 50% of the length of WINDOW if not specified.
%
%  [L,P,F] = TOD(X,WINDOW,WEIGHT,NOVERLAP,DT) uses the time step DT
%  between consecutive snapshots to determine a physical frequency F.
%
%  [L,P,F] = TOD(X,WINDOW,WEIGHT,NOVERLAP,DT,OPTS) specifies options:
%  OPTS.Q           function that takes inputs (q1,q2) and returns the
%                   nonlinear term of the equation [function |
%                   {@(q1,q2) permute(q1.*q2,[2 1 3])}]
%  OPTS.LHS         function that takes input (q) and returns the left hand
%                   side of the equation [function | {@(q) q}]
%  OPTS.nmode       number of modes (ranks) per triad to store [integer |
%                   {nBlks}]
%  OPTS.threads     number of threads if MATLAB supports parallel
%                   functionality; if not, or if OPTS.threads<=1, the
%                   code runs serially [integer | {maxNumCompThreads}]
%  OPTS.isreal      reality of data, which determines regions of the
%                   bispectrum to compute, see figure caption above
%                   [logical | {isreal(X)}]
%  OPTS.precision   compute TOD in single precision ['single' | 'double' |
%                   {class(X)}]
%  OPTS.mean        provide a mean that is subtracted from each
%                   snapshot [array of size X | {0}]
%  OPTS.nfreq       restrict computation to |l|,|k|,|n|<=OPTS.nfreq
%                   [integer | {all}] 
%
%  [L,P,F,T] = TOD(...) returns the modal energy budget in T. T has the
%  same dimensions as L.
%
%  [L,P,F,T,A] = TOD(...) returns the expansion coefficients in A. The
%  matrices of coefficients associated with the convective and recipient
%  modes are stored respectively in A(1,l,n,:,:) and A(2,l,n,:,:).
%
%  [L,P,F,T,A,Xi] = TOD(...) returns the donor modes in Xi{l,n}(1,...) and
%  catalyst modes in Xi{l,n}(2,...), where l and n are the same indices that
%  index into L, P, T, and A.
%
%  References:
%   [1] Yeung, B., Chu, T., and Schmidt, O. T., Triadic orthogonal
%       decomposition reveals nonlinearity in fluid flows, J. Fluid Mech. 2026
%       DOI 10.1017/jfm.2026.11183
%   [2] Schmidt, O. T., Bispectral mode decomposition of nonlinear flows, 
%       Nonlinear Dyn. 2020
%       DOI 10.1007/s11071-020-06037-z
%
% B. Yeung (byeung@ucsd.edu), T. Chu, and O. T. Schmidt
% Original script: BY and TC, based on BMD code [2] by OTS
% Last revision:   16-Apr-2026 (BY)

if nargin==6
    opts = varargin{5};
    if ~isfield(opts,'isreal')
        opts.isreal = isreal(X);
    end
    if isfield(opts,'precision')
        if strncmpi(opts.precision,'single',6)
            opts.precision = 'single';
        elseif strncmpi(opts.precision,'double',6)
            opts.precision = 'double';
        else
            opts.precision = class(X);
        end
    else
        opts.precision = class(X);
    end
else
    opts.isreal = isreal(X);
    opts.precision = class(X);
end

% get problem dimensions
dim     = size(X);
nt      = dim(1);
nVar    = dim(2);
nx      = prod(dim(3:end));

% get default spectral estimation parameters and options
[window,weight,nOvlp,dt,nDFT,nBlks] = parser(nt,nx,varargin{:});

% determine correction for FFT window gain
winWeight   = 1/mean(window);

% nonlinear quadratic nonlinearity
if isfield(opts,'Q')
    Q_nonlinear = opts.Q ;
else
    Q_nonlinear = @(q1,q2) permute(q1.*q2,[2 1 3]);
end

% lhs
if isfield(opts,'LHS')
    LHS = opts.LHS;
else
    LHS = @(q) q;
end

% use long-time mean if provided
blk_mean = false;
if isfield(opts,'mean')
    if strcmp(opts.mean,'none')
    X_mean      = zeros(nVar,1);
    mean_name   = 'zero';
    elseif strcmp(opts.mean,'blockwise')
    blk_mean = true;
    X_mean      = zeros(nVar,1);
    mean_name   = 'blockwise mean';
    else
    X_mean      = opts.mean(1:nVar,:);
    mean_name   = 'provided long-time mean';    
    end
else
    X_mean      = zeros(nVar,1);
    mean_name   = 'zero'; 
end

disp(['Mean                      : ' mean_name]);

% obtain frequency axis
[f,nFreq,includeTriad,f_idx,fk_idx,fl_idx,fn_idx] = faxes(nDFT,dt,opts);
nTriads = nnz(includeTriad);

% loop over number of blocks and generate Fourier realizations
disp(' ')
disp('Calculating temporal DFT')
disp('------------------------------------')
Q_hat = zeros(nFreq,nVar,nx,nBlks,opts.precision);
for iBlk = 1:nBlks
    % get time index for present block
    offset                  = min((iBlk-1)*(nDFT-nOvlp)+nDFT,nt)-nDFT;
    timeIdx                 = (1:nDFT) + offset;
    disp(['block ' num2str(iBlk) '/' num2str(nBlks) ' (' ...
        num2str(timeIdx(1)) ':' num2str(timeIdx(end)) ')'])
    
    for iVar = 1:nVar
        Q_blk          = bsxfun(@minus,squeeze(X(timeIdx,iVar,:)),squeeze(X_mean(iVar,:)));
        if blk_mean
            Q_blk_tmp = Q_blk;
            Q_blk = Q_blk-mean(Q_blk,1);
        end
        Q_blk          = bsxfun(@times,Q_blk,window);
        Q_blk_hat      = winWeight/nDFT*fft(Q_blk);
        if blk_mean, Q_blk_hat(1,:,:) = winWeight/nDFT*sum(Q_blk_tmp.*window,1); end
        Q_blk_hat      = fftshift(Q_blk_hat,1);
        Q_hat(:,iVar,:,iBlk) = Q_blk_hat;
    end
end
clear X Q_blk Q_blk_hat X_mean

nState = size(LHS(squeeze(Q_hat(1,:,:,:))),1);
weights = repmat(weight,nState,1);

% preallocate numeric and cell arrays for bispectra and modes respectively
if ~isfield(opts,'nmode'), opts.nmode = nBlks; end
L = nan(nFreq,nFreq,nBlks,opts.precision);
nout = nargout;
if nout>3
    T = nan(nFreq,nFreq,nBlks,opts.precision);
end
P = cell(nFreq,nFreq);
if nout>4, A = nan(2,nFreq,nFreq,nBlks,nBlks,opts.precision); end
if nout>5, Xi = cell(nFreq,nFreq); end

% loop over all triads and calculate TOD
disp(' ')
disp('Calculating TOD')
disp('------------------------------------')

if ~isfield(opts,'threads'), opts.threads = maxNumCompThreads; end
if supportsParallel && opts.threads>1
    % if a pool already exists, close it
    if ~isempty(gcp('nocreate')), delete(gcp('nocreate')); end
    parpool('threads',opts.threads);

    % calculate TOD in parallel mode
    doParallelTOD

    delete(gcp('nocreate'))
else
    % calculate TOD in serial mode
    doSerialTOD
end

P   = cellfun(@(p) reshapecell(p,[2 dim(3:end) nState opts.nmode]),P,'uniformoutput',false);
if nout>5, Xi = cellfun(@(p) reshapecell(p,[2 dim(3:end) nState opts.nmode]),Xi,'uniformoutput',false); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nested functions
function doParallelTOD
spmd
    triadIndices = getLocalIndices(nTriads);
for i = triadIndices
    disp(['(' num2str(f_idx(fk_idx(i))) ',' num2str(f_idx(fl_idx(i))) ',' num2str(f_idx(fn_idx(i))) ') (' num2str(i) '/' num2str(nTriads) ')'])
    Q_hat_n = reshape(permute(LHS(squeeze(Q_hat(fn_idx(i),:,:,:))),[2 1 3]),nState*nx,nBlks);
    Q_hat_kl = reshape(Q_nonlinear(squeeze(Q_hat(fk_idx(i),:,:,:)),squeeze(Q_hat(fl_idx(i),:,:,:))),nState*nx,nBlks);
    Q_hat_l = reshape(permute(LHS(squeeze(Q_hat(fl_idx(i),:,:,:))),[2 1 3]),nState*nx,nBlks);
    Q_hat_k = reshape(permute(LHS(squeeze(Q_hat(fk_idx(i),:,:,:))),[2 1 3]),nState*nx,nBlks);
        
    [U,s,V] = todAlgorithm(Q_hat_n,Q_hat_kl,weights,nBlks);
    u = U(:,1:opts.nmode);
    v = V(:,1:opts.nmode);

    P{fl_idx(i),fn_idx(i)}(1,:,:) = u; % convective modes (l->n)
    P{fl_idx(i),fn_idx(i)}(2,:,:) = v; % recipient modes (n)
    L(fl_idx(i),fn_idx(i),:) = s; % mode bispectra (singular values)

    if nout>3
        % integral modal energy transfers
        T(fl_idx(i),fn_idx(i),:) = s.'.*dot(V,bsxfun(@times,U,weights));
        
        if nout>4
            % expansion coefficients
            A(1,fl_idx(i),fn_idx(i),:,:) = U'*(Q_hat_kl.*weights);
            A(2,fl_idx(i),fn_idx(i),:,:) = V'*(Q_hat_n.*weights);

            % donor and catalyst modes
            if nout>5
                A_this = squeeze(A(2,fl_idx(i),fn_idx(i),:,:));
                % donor (l), catalyst (k=n-l)
                donor_mode = Q_hat_l*A_this'*diag(1./s)/nBlks;
                catalyst_mode = Q_hat_k*A_this'*diag(1./s)/nBlks;
                Xi{fl_idx(i),fn_idx(i)}(1,:,:) = donor_mode(:,1:opts.nmode);
                Xi{fl_idx(i),fn_idx(i)}(2,:,:) = catalyst_mode(:,1:opts.nmode);
            end
        end
    end
end
L = spmdReduce(@addnan,L,1);
P = spmdReduce(@addcell,P,1);
if nout>3
    T = spmdReduce(@addnan,T,1);
    if nout>4
        A = spmdReduce(@addnan,A,1);
        if nout>5
            Xi = spmdReduce(@addcell,Xi,1);
        end
    end
end
end
L = gather(L);
P = gather(P);
L = L{1};
P = P{1};
if nout>3
    T = gather(T);
    T = T{1};
    if nout>4
        A = gather(A);
        A = A{1};
        if nout>5
            Xi = gather(Xi);
            Xi = Xi{1};
        end
    end
end
end

function doSerialTOD
for i = 1:nTriads
    disp(['(' num2str(f_idx(fk_idx(i))) ',' num2str(f_idx(fl_idx(i))) ',' num2str(f_idx(fn_idx(i))) ') (' num2str(i) '/' num2str(nTriads) ')'])
    Q_hat_n = reshape(permute(LHS(squeeze(Q_hat(fn_idx(i),:,:,:))),[2 1 3]),nState*nx,nBlks);
    Q_hat_kl = reshape(Q_nonlinear(squeeze(Q_hat(fk_idx(i),:,:,:)),squeeze(Q_hat(fl_idx(i),:,:,:))),nState*nx,nBlks);
    Q_hat_l = reshape(permute(LHS(squeeze(Q_hat(fl_idx(i),:,:,:))),[2 1 3]),nState*nx,nBlks);
    Q_hat_k = reshape(permute(LHS(squeeze(Q_hat(fk_idx(i),:,:,:))),[2 1 3]),nState*nx,nBlks);
        
    [U,s,V] = todAlgorithm(Q_hat_n,Q_hat_kl,weights,nBlks);
    u = U(:,1:opts.nmode);
    v = V(:,1:opts.nmode);

    P{fl_idx(i),fn_idx(i)}(1,:,:) = u; % convective modes (l->n)
    P{fl_idx(i),fn_idx(i)}(2,:,:) = v; % recipient modes (n)
    L(fl_idx(i),fn_idx(i),:) = s; % mode bispectra (singular values)

    if nout>3
        % integral modal energy transfers
        T(fl_idx(i),fn_idx(i),:) = s.'.*dot(V,bsxfun(@times,U,weights));
        
        if nout>4
            % expansion coefficients
            A(1,fl_idx(i),fn_idx(i),:,:) = U'*(Q_hat_kl.*weights);
            A(2,fl_idx(i),fn_idx(i),:,:) = V'*(Q_hat_n.*weights);

            % donor and catalyst modes
            if nout>5
                A_this = squeeze(A(2,fl_idx(i),fn_idx(i),:,:));
                % donor (l), catalyst (k=n-l)
                donor_mode = Q_hat_l*A_this'*diag(1./s)/nBlks;
                catalyst_mode = Q_hat_k*A_this'*diag(1./s)/nBlks;
                Xi{fl_idx(i),fn_idx(i)}(1,:,:) = donor_mode(:,1:opts.nmode);
                Xi{fl_idx(i),fn_idx(i)}(2,:,:) = catalyst_mode(:,1:opts.nmode);
            end
        end
    end
end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [window,weight,nOvlp,dt,nDFT,nBlks] = parser(nt,nx,varargin)
% PARSER Parser for parameters

% read input arguments from cell array
window = []; weight = []; nOvlp = []; dt = [];
nvarargin = length(varargin);

if nvarargin >= 1
    window = varargin{1};
    if nvarargin >= 2
        weight   = varargin{2};
        if nvarargin >= 3
            nOvlp   = varargin{3};
            if nvarargin >= 4
                dt      = varargin{4};
            end
        end
    end
end

window = window(:); weight = weight(:);

% check arguments and determine default spectral estimation parameters
% window size and type
if isempty(window)
    nDFT        = 2^floor(log2(nt/5)); if nDFT>256, nDFT=256; end
    window      = hammwin(nDFT);
    window_name = 'Hamming';
elseif length(window)==1
    nDFT        = window;
    window      = hammwin(window);
    window_name = 'Hamming';
elseif length(window) == 2^nextpow2(length(window))
    nDFT        = length(window);
    window_name = 'user specified';
else
    nDFT        = length(window);
    window_name = 'user specified';
end

% block overlap
if isempty(nOvlp)
    nOvlp = floor(nDFT/2);
elseif nOvlp > nDFT-1
    error('Overlap too large.')
end

% time step between consecutive snapshots
if isempty(dt)
    dt = 1/nDFT;
end

% inner product weight
if isempty(weight)
    weight      = ones(nx,1);
    weight_name = 'uniform';
elseif numel(weight) ~= nx
    error('Weights must have the same spatial dimensions as data.');
else
    weight_name = 'user specified';
end

% number of blocks
nBlks = floor((nt-nOvlp)/(nDFT-nOvlp));

% test feasibility
if nDFT < 4 || nBlks < 2
    error('Spectral estimation parameters not meaningful.');
end

% display parameter summary
disp(' ')
disp('TOD parameters')
disp('------------------------------------')
disp(['No. of snaphots per block : ' num2str(nDFT)])
disp(['Block overlap             : ' num2str(nOvlp)])
disp(['No. of blocks             : ' num2str(nBlks)])
disp(['Windowing fct. (time)     : ' window_name])
disp(['Weighting fct. (space)    : ' weight_name])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [window] = hammwin(N)
% HAMMWIN standard Hamming window of lenght N
window = 0.54-0.46*cos(2*pi*(0:N-1)/(N-1))';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = supportsParallel
% brute-force check of parallel support
try
    gcp('nocreate');
    tf = true;
catch
    tf = false;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,s,V] = todAlgorithm(Q_hat_n,Q_hat_kl,weights,nBlks)
[U,S,V] = lowrankSVD(bsxfun(@times,Q_hat_n,sqrt(weights))'/nBlks,bsxfun(@times,Q_hat_kl,sqrt(weights)));
U = bsxfun(@times,U,1./sqrt(weights));
V = bsxfun(@times,V,1./sqrt(weights));
s = diag(S);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,S,V] = lowrankSVD(X,Q3)
[Q,R] = qr(Q3,0);
[U,S,V] = sirovichSVD(R*X);
U = Q*U;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,S,V] = sirovichSVD(X)
[U,L] = eig(X*X','vector');
[L,idx] = sort(L,'descend');
U = U(:,idx);
sqrtL = sqrt(L);
S = diag(sqrtL);
V = X'*U*diag(1./sqrtL);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,nFreq,includeTriad,f_idx,fk_idx,fl_idx,fn_idx] = faxes(nDFT,dt,opts)
% FAXES obtain frequency axes and indices
f_idx = (0:nDFT-1);
if mod(nDFT,2)==0
    f_idx(nDFT/2+1:end)     = f_idx(nDFT/2+1:end)-nDFT;
else
    f_idx((nDFT+1)/2+1:end) = f_idx((nDFT+1)/2+1:end)-nDFT;
end
f_idx   = fftshift(f_idx);
f       = f_idx/dt/nDFT;
nFreq   = numel(f_idx);
f_idx_min   = f_idx(1);
f_idx_max   = f_idx(end);
if isfield(opts,'nfreq')
    if opts.nfreq<f_idx_max
        f_idx_max = opts.nfreq;
        f_idx_min = -f_idx_max;
    end
end

includeTriad = ones(nFreq,nFreq,'logical');
f0_idx = floor(nFreq/2)+1; % index of f=0
[fl_idx,fn_idx] = ndgrid(1:nFreq);
fk_idx = toeplitz(f0_idx:-1:f0_idx-nFreq+1,f0_idx:f0_idx+nFreq-1);

fl = f_idx(fl_idx);
fn = f_idx(fn_idx);
fk = fn-fl;

% if real data, principal region is top half of (fl,fn) plane
if opts.isreal
    includeTriad(fn<0) = 0;
end

% truncate bispectrum when fk, fl, or fn exceeds Nyquist or user-specified
% max frequency
includeTriad(fk<f_idx_min | fk>f_idx_max | fl<f_idx_min | fl>f_idx_max | fn<f_idx_min | fn>f_idx_max) = 0;

% truncate frequency indices
fk_idx = fk_idx(includeTriad);
fl_idx = fl_idx(includeTriad);
fn_idx = fn_idx(includeTriad);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = addnan(a,b)
%ADDNAN consolidates non-NaN elements in arrays a and b into b

b(~isnan(a)) = a(~isnan(a));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = addcell(a,b)
%ADDCELL consolidates non-empty cells in cell arrays a and b into b

b(~cellfun(@isempty,a)) = a(~cellfun(@isempty,a));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = reshapecell(a,output_dim)
%RESHAPECELL reshapes non-empty a to dimensions output_dim

if ~isempty(a)
    a = reshape(a,output_dim);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localTriadIndices = getLocalIndices(nTriads)
%LOCALTRIADINDICES returns triad indices local to each thread
nTriads_local = floor(nTriads/numlabs);
if labindex~=numlabs
    localTriadIndices = (labindex-1)*nTriads_local+(1:nTriads_local);
else
    localTriadIndices = ((labindex-1)*nTriads_local+1):nTriads;
end
end