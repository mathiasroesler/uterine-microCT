function STPipeline(data_folder, diffusion, structure_tensor, streamlines)
%STPipeline Runs the structure tensor analysis pipeline.
%
%   Input:
%    - data_folder, name of the data folder to use located in the
%    microCT/data folder.
%    - diffusion, true if the DiffusionTissueExtrapolation code should be
%    run, default value is true.
%    - structure_tensor, true if the StructureTensorExtraction code should
%    be run, default value is true.
%    - streamlines, true if the ComputeStreamlines code should be run,
%    default value is true.
%
%   Return:
%
if nargin < 4
    streamlines = true; 
end

if nargin < 3
    structure_tensor = true;
end

if nargin < 2
    diffusion = true;
end

%% General parameters
base_dir = join([getenv("HOME"), "Documents/phd/microCT/data"], '/');
src_dir = join([base_dir, data_folder, "downsampled/ST"], '/');

[file_template, extension] = loadParams(src_dir + '/ST.params');

if diffusion
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % DiffusionTissueExtrapolation.m
    %
    % This script loads in a masked tissue image stack and identifies the
    % tissue as having values > threshold. The tissue mask is heavily
    % eroded to ensure that tissue intensity values provide boundary
    % conditions for diffusion extrapolation into non-tissue regions.
    %
    % The diffusion extrapolation iterates until the mean of the differences
    % in intensity between voxels a specified distance from the tissue
    % and the nearest tissue values crosses zero.
    %
    % Written by: Mark Trew, November 2011.
    % Modified by: Mark Trew, March 2021
    % Modified by: Mark Trew, September 2021
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % User defined parameters
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    img_input_dir = src_dir + '/masked/';
    img_output_dir = src_dir + '/extrapolated/';

    % Tissue mask erosion threshold and radius
    TissueMaskErosionThreshold = 10;
    TissueMaskErosionRadius = 2;

    % Set the tissue boundary diffusion testing distance (voxels) - points
    % at this distance from the tissue will be used to assess the stopping
    % criteria for the diffusion iterations.
    TissueBoundaryDiffusionTestingDistance = 2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Set up and load in image stack
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('... loading images ...\n');
    img_paths = getImagePaths(img_input_dir, extension);
    I = loadImageStack(img_paths);
    [Nj, Ni, Nk] = size(I);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Identify tissue mask and erode.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('...mask and erosion...\n');
    % 1. Identify tissue mask
    IB = (I >= TissueMaskErosionThreshold);
    % 2. erode tissue mask - the erosion filter is spherical with diameter
    %    2*R+1 voxels
    R = TissueMaskErosionRadius;
    NHood = zeros(2*R+1,2*R+1,2*R+1);
    NHood(R+1,R+1,R+1) = 1;
    D = bwdist(NHood);
    NHood = (D <= R);
    IB = imerode(IB,NHood);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Find voxels a specified distance from the tissue. Store nearest tissue
    % values.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('...finding Euclidian distances...\n');
    DT = TissueBoundaryDiffusionTestingDistance;
    [D,L]=bwdist(IB);
    TIdx = find(D >= DT-0.5 & D <= DT+0.5);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Set up and perform diffusion.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set timing initialization
    trun0 = cputime;

    % Identify the points that will be used as boundary conditions
    % on the diffusion
    FixedPoints = IB;

    % 3D 7 point diffusion filter
    F = zeros(3,3,3);
    F(2,2,1)=1;
    F(2,1,2)=1; F(1,2,2)=1; F(3,2,2)=1; F(2,3,2)=1; F(2,2,2)=-6;
    F(2,2,3)=1;

    % Set up a 3D 7 point diffusion filter and then construct a circulant
    % matrix.

    % Initial conditions on diffused image DI. Delta is the effective
    % "time step". If it is too large the diffusion will be unstable.
    % Clear the I and IB arrays as they are no longer needed.
    DI = I.*cast(IB,'uint8'); clear I  IB  D;
    Delta = 0.3;

    % Iterate over diffusion steps while the sign of the mean difference
    % between the diffused image intensity DT voxels from the tissue and the
    % nearest boundary condition is less than zero. The L array is the index of
    % the original closest non-zero pixel to the test group, i.e. the nearest
    % boundary condition.
    % Note: this algorithm could be rewritten to potentially be faster.
    % Using the same tools as for structure tensor analysis the F could be
    % set up as a circulant matrix. Then the update operation would be:
    % DI = Delta*F*DI + DI = (Delta*F+I)*DI = S*DI
    % Here S becomes a new circulant matrix and the S*DI operation simply
    % becomes a convolution operation and should be relatively fast.
    n=1;
    md = -1;
    while  md <= 0
        md = mean(double(DI(TIdx))-double(DI(L(TIdx))));
        if ~mod(n,10), fprintf('... n: %d, mean diff: %f\n',n,md); end
        % apply diffusion filter to current state of diffused image. DI is
        % exterior padded with replicated values to enforce zero flux b.c.
        DIhat = imfilter(padarray(DI,[1,1,1],'replicate'),F);
        DIhat = DIhat(2:Nj+1,2:Ni+1,2:Nk+1);
        % apply fixed boundary conditions
        DIhat(FixedPoints) = 0;
        % update the diffused image
        DI = Delta*DIhat + DI;
        % update counter
        n = n+1;
    end

    trun1 = cputime;
    fprintf(' Total time for %d its: %f\n',n,trun1-trun0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Write boundary smoothed image files.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    saveImageStack(DI, img_output_dir, file_template);

end

if structure_tensor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % StructureTensorExtraction.m
    %
    % This script loads in a stack of images (with isotropic voxels) and
    % computes the components of the structure tensor (the outer product of
    % the intensity gradient vectors). The components are smoothed to a
    % sequence of 1/2 resolutions and these are written to file.
    %
    % Modified by: Mark Trew, September 2021.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Timing
    ttotalprocess0 = cputime;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % User defined parameters
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    img_input_dir = src_dir + '/extrapolated/';
    DataOutput = src_dir + '/binary/';

    % Set the derivative and smoothing template voxel widths
    DerivativeTemplateWidth = 3;
    SmoothingTemplateWidth = 3;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Set up and load in image stack
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Load in image set
    fprintf('... loading images ...\n');
    img_paths = getImagePaths(img_input_dir, extension);
    I = loadImageStack(img_paths);
    [Nj, Ni, Nk] = size(I);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Rearrange and pad the image stack
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Permute the array to make the x/i index the most rapidly varying.
    % This process is necessary to ensure the i,j,k indices of the matlab
    % array correspond with the xyz coordinates of the texture orientation
    % vectors.
    fprintf('... permuting array ...\n');
    I = permute(I,[2,1,3]);

    % Pad the image around each edge - reflect image in padding
    fprintf('... padding array ...\n');
    IPad = padarray(I,[2,2,2],'symmetric'); clear I;

    % Reshape padded permuted equalized image to a 1D array
    IPad = reshape(IPad,(Ni+4)*(Nj+4)*(Nk+4),1);

    IPad = single(IPad)/single(max(IPad));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Derivative filters
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up the 1D gradient derivative weight filters
    fprintf('... computing gradient weight templates ...\n');
    [Wx,Wy,Wz] = ConstructDerivativeFilters(Ni,Nj,Nk,DerivativeTemplateWidth);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % FFT calculations
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute fft of Image
    fprintf('... computing FFT of image ...\n');
    t10 = clock;
    IPadf = single(fft(IPad)); %clear IPad;
    t1 = etime(clock,t10); fprintf(' fft time for IPad: %0.2f sec\n',t1);

    % Compute fft of 1st derivative weights
    fprintf('... computing FFTs of 1st derivative weights ...\n');
    fprintf('... performing FFT convolution for 1st derivatives ...\n');
    t10 = clock;
    Wxf = single(fft(Wx)); clear Wx;
    Di = double(ifft(IPadf .* Wxf)); clear Wxf;% FFT convolution
    t1 = etime(clock,t10); fprintf(' FFT Dx Convolution time for Wx: %0.2f sec\n',t1);
    t10 = clock;
    Wyf = single(fft(Wy)); clear Wy;
    Dj = double(ifft(IPadf .* Wyf)); clear Wyf;% FFT convolution
    t1 = etime(clock,t10); fprintf(' FFT Dy Convolution time for Wy: %0.2f sec\n',t1);
    t10 = clock;
    Wzf = single(fft(Wz)); clear Wz;
    Dk = double(ifft(IPadf .* Wzf)); clear Wzf;% FFT convolution
    t1 = etime(clock,t10); fprintf(' FFT Dz Convolution time for Wz: %0.2f sec\n',t1);

    clear IPadf;

    % Useful data ranges
    PL = 2; % valid derivative padding level
    SBi = 1+1*PL:Ni+4-1*PL;
    SBj = 1+1*PL:Nj+4-1*PL;
    SBk = 1+1*PL:Nk+4-1*PL;

    % Raw index locations
    [SI0,SJ0,SK0] = ndgrid(1:length(SBi),1:length(SBj),1:length(SBk));

    % Report size of zero level data
    fprintf('... Zero level data dimensions: (%d,%d,%d)\n',size(SI0));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Perform smoothing
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Set up for first level smoothing
    t0 = clock;
    fprintf('... First level ST construction and smoothing ...\n');

    fprintf('     - Jii ...\n');
    Jii = Di.*Di; Jii = reshape(Jii,[(Ni+4),(Nj+4),(Nk+4)]);
    [Sii1,SI1,SJ1,SK1] = MultigridAveraging(Jii(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);

    fprintf('     - Jij ...\n');
    Jij = Di.*Dj; Jij = reshape(Jij,[(Ni+4),(Nj+4),(Nk+4)]);
    [Sij1,SI1,SJ1,SK1] = MultigridAveraging(Jij(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);

    fprintf('     - Jik ...\n');
    Jik = Di.*Dk; Jik = reshape(Jik,[(Ni+4),(Nj+4),(Nk+4)]); %clear Di;
    [Sik1,SI1,SJ1,SK1] = MultigridAveraging(Jik(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);

    fprintf('     - Jjj ...\n');
    Jjj = Dj.*Dj; Jjj = reshape(Jjj,[(Ni+4),(Nj+4),(Nk+4)]);
    [Sjj1,SI1,SJ1,SK1] = MultigridAveraging(Jjj(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);

    fprintf('     - Jjk ...\n');
    Jjk = Dj.*Dk; Jjk = reshape(Jjk,[(Ni+4),(Nj+4),(Nk+4)]); %clear Dj;
    [Sjk1,SI1,SJ1,SK1] = MultigridAveraging(Jjk(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);

    fprintf('     - Jkk ...\n');
    Jkk = Dk.*Dk; Jkk = reshape(Jkk,[(Ni+4),(Nj+4),(Nk+4)]); %clear Dk;
    [Skk1,SI1,SJ1,SK1] = MultigridAveraging(Jkk(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);

    t1 = etime(clock,t0); fprintf(' First level ST construction and smoothing time: %0.2f sec\n',t1);

    fprintf('... First level data dimensions: (%d,%d,%d)\n',size(SI1));

    % Smooth to second level using multigrid binomial averaging
    fprintf('... Second level smoothing ...\n');
    t0 = clock;
    [Sii2,SI2,SJ2,SK2] = MultigridAveraging(Sii1,SI1,SJ1,SK1,SmoothingTemplateWidth);
    [Sij2,SI2,SJ2,SK2] = MultigridAveraging(Sij1,SI1,SJ1,SK1,SmoothingTemplateWidth);
    [Sik2,SI2,SJ2,SK2] = MultigridAveraging(Sik1,SI1,SJ1,SK1,SmoothingTemplateWidth);
    [Sjj2,SI2,SJ2,SK2] = MultigridAveraging(Sjj1,SI1,SJ1,SK1,SmoothingTemplateWidth);
    [Sjk2,SI2,SJ2,SK2] = MultigridAveraging(Sjk1,SI1,SJ1,SK1,SmoothingTemplateWidth);
    [Skk2,SI2,SJ2,SK2] = MultigridAveraging(Skk1,SI1,SJ1,SK1,SmoothingTemplateWidth);
    t1 = etime(clock,t0); fprintf(' Second level smoothing time: %0.2f sec\n',t1);
    fprintf('... Second level data dimensions: (%d,%d,%d)\n',size(SI2));

    % Smooth to third level using multigrid binomial averaging
    fprintf('... Third level smoothing ...\n');
    t0 = clock;
    [Sii3,SI3,SJ3,SK3] = MultigridAveraging(Sii2,SI2,SJ2,SK2,SmoothingTemplateWidth);
    [Sij3,SI3,SJ3,SK3] = MultigridAveraging(Sij2,SI2,SJ2,SK2,SmoothingTemplateWidth);
    [Sik3,SI3,SJ3,SK3] = MultigridAveraging(Sik2,SI2,SJ2,SK2,SmoothingTemplateWidth);
    [Sjj3,SI3,SJ3,SK3] = MultigridAveraging(Sjj2,SI2,SJ2,SK2,SmoothingTemplateWidth);
    [Sjk3,SI3,SJ3,SK3] = MultigridAveraging(Sjk2,SI2,SJ2,SK2,SmoothingTemplateWidth);
    [Skk3,SI3,SJ3,SK3] = MultigridAveraging(Skk2,SI2,SJ2,SK2,SmoothingTemplateWidth);
    t1 = etime(clock,t0); fprintf(' Third level smoothing time: %0.2f sec\n',t1);
    fprintf('... Third level data dimensions: (%d,%d,%d)\n',size(SI3));

    % Smooth to fourth level using multigrid binomial averaging
    fprintf('... Fourth level smoothing ...\n');
    t0 = clock;
    [Sii4,SI4,SJ4,SK4] = MultigridAveraging(Sii3,SI3,SJ3,SK3,SmoothingTemplateWidth);
    [Sij4,SI4,SJ4,SK4] = MultigridAveraging(Sij3,SI3,SJ3,SK3,SmoothingTemplateWidth);
    [Sik4,SI4,SJ4,SK4] = MultigridAveraging(Sik3,SI3,SJ3,SK3,SmoothingTemplateWidth);
    [Sjj4,SI4,SJ4,SK4] = MultigridAveraging(Sjj3,SI3,SJ3,SK3,SmoothingTemplateWidth);
    [Sjk4,SI4,SJ4,SK4] = MultigridAveraging(Sjk3,SI3,SJ3,SK3,SmoothingTemplateWidth);
    [Skk4,SI4,SJ4,SK4] = MultigridAveraging(Skk3,SI3,SJ3,SK3,SmoothingTemplateWidth);
    t1 = etime(clock,t0); fprintf(' Fourth level smoothing time: %0.2f sec\n',t1);
    fprintf('... Fourth level data dimensions: (%d,%d,%d)\n',size(SI4));

    % Smooth to fifth level using multigrid binomial averaging
    fprintf('... Fifth level smoothing ...\n');
    t0 = clock;
    [Sii5,SI5,SJ5,SK5] = MultigridAveraging(Sii4,SI4,SJ4,SK4,SmoothingTemplateWidth);
    [Sij5,SI5,SJ5,SK5] = MultigridAveraging(Sij4,SI4,SJ4,SK4,SmoothingTemplateWidth);
    [Sik5,SI5,SJ5,SK5] = MultigridAveraging(Sik4,SI4,SJ4,SK4,SmoothingTemplateWidth);
    [Sjj5,SI5,SJ5,SK5] = MultigridAveraging(Sjj4,SI4,SJ4,SK4,SmoothingTemplateWidth);
    [Sjk5,SI5,SJ5,SK5] = MultigridAveraging(Sjk4,SI4,SJ4,SK4,SmoothingTemplateWidth);
    [Skk5,SI5,SJ5,SK5] = MultigridAveraging(Skk4,SI4,SJ4,SK4,SmoothingTemplateWidth);
    t1 = etime(clock,t0); fprintf(' Fifth level smoothing time: %0.2f sec\n',t1);
    fprintf('... Fifth level data dimensions: (%d,%d,%d)\n',size(SI5));

    % Smooth to sixth level using multigrid binomial averaging
    fprintf('... Sixth level smoothing ...\n');
    t0 = clock;
    [Sii6,SI6,SJ6,SK6] = MultigridAveraging(Sii5,SI5,SJ5,SK5,SmoothingTemplateWidth);
    [Sij6,SI6,SJ6,SK6] = MultigridAveraging(Sij5,SI5,SJ5,SK5,SmoothingTemplateWidth);
    [Sik6,SI6,SJ6,SK6] = MultigridAveraging(Sik5,SI5,SJ5,SK5,SmoothingTemplateWidth);
    [Sjj6,SI6,SJ6,SK6] = MultigridAveraging(Sjj5,SI5,SJ5,SK5,SmoothingTemplateWidth);
    [Sjk6,SI6,SJ6,SK6] = MultigridAveraging(Sjk5,SI5,SJ5,SK5,SmoothingTemplateWidth);
    [Skk6,SI6,SJ6,SK6] = MultigridAveraging(Skk5,SI5,SJ5,SK5,SmoothingTemplateWidth);

    % Smooth to seventh level using multigrid binomial averaging
    % fprintf('... Seventh level smoothing ...\n');
    % t0 = clock;
    % [Sii7,SI7,SJ7,SK7] = MultigridAveraging(Sii6,SI6,SJ6,SK6,SmoothingTemplateWidth);
    % [Sij7,SI7,SJ7,SK7] = MultigridAveraging(Sij6,SI6,SJ6,SK6,SmoothingTemplateWidth);
    % [Sik7,SI7,SJ7,SK7] = MultigridAveraging(Sik6,SI6,SJ6,SK6,SmoothingTemplateWidth);
    % [Sjj7,SI7,SJ7,SK7] = MultigridAveraging(Sjj6,SI6,SJ6,SK6,SmoothingTemplateWidth);
    % [Sjk7,SI7,SJ7,SK7] = MultigridAveraging(Sjk6,SI6,SJ6,SK6,SmoothingTemplateWidth);
    % [Skk7,SI7,SJ7,SK7] = MultigridAveraging(Skk6,SI6,SJ6,SK6,SmoothingTemplateWidth);
    % t1 = etime(clock,t0); fprintf(' Seventh level smoothing time: %0.2f sec\n',t1);
    % fprintf('... Seventh level data dimensions: (%d,%d,%d)\n',size(SI7));

    ttotalprocess1 = cputime; fprintf(' *** total processing time: %0.2f sec\n',ttotalprocess1-ttotalprocess0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Write binary data files
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('... Writing binary data files ...\n');

    fid = fopen(sprintf('%sS2.bin',DataOutput),'wb');
    fwrite(fid,[size(Sii2,1),size(Sii2,2),size(Sii2,3)],'uint16');
    fwrite(fid,reshape(Sii2,numel(Sii2),1),'double');
    fwrite(fid,reshape(Sij2,numel(Sij2),1),'double');
    fwrite(fid,reshape(Sik2,numel(Sik2),1),'double');
    fwrite(fid,reshape(Sjj2,numel(Sjj2),1),'double');
    fwrite(fid,reshape(Sjk2,numel(Sjk2),1),'double');
    fwrite(fid,reshape(Skk2,numel(Skk2),1),'double');
    fclose(fid);

    fid = fopen(sprintf('%sS3.bin',DataOutput),'wb');
    fwrite(fid,[size(Sii3,1),size(Sii3,2),size(Sii3,3)],'uint16');
    fwrite(fid,reshape(Sii3,numel(Sii3),1),'double');
    fwrite(fid,reshape(Sij3,numel(Sij3),1),'double');
    fwrite(fid,reshape(Sik3,numel(Sik3),1),'double');
    fwrite(fid,reshape(Sjj3,numel(Sjj3),1),'double');
    fwrite(fid,reshape(Sjk3,numel(Sjk3),1),'double');
    fwrite(fid,reshape(Skk3,numel(Skk3),1),'double');
    fclose(fid);

    fid = fopen(sprintf('%sS4.bin',DataOutput),'wb');
    fwrite(fid,[size(Sii4,1),size(Sii4,2),size(Sii4,3)],'uint16');
    fwrite(fid,reshape(Sii4,numel(Sii4),1),'double');
    fwrite(fid,reshape(Sij4,numel(Sij4),1),'double');
    fwrite(fid,reshape(Sik4,numel(Sik4),1),'double');
    fwrite(fid,reshape(Sjj4,numel(Sjj4),1),'double');
    fwrite(fid,reshape(Sjk4,numel(Sjk4),1),'double');
    fwrite(fid,reshape(Skk4,numel(Skk4),1),'double');
    fclose(fid);

    fid = fopen(sprintf('%sS5.bin',DataOutput),'wb');
    fwrite(fid,[size(Sii5,1),size(Sii5,2),size(Sii5,3)],'uint16');
    fwrite(fid,reshape(Sii5,numel(Sii5),1),'double');
    fwrite(fid,reshape(Sij5,numel(Sij5),1),'double');
    fwrite(fid,reshape(Sik5,numel(Sik5),1),'double');
    fwrite(fid,reshape(Sjj5,numel(Sjj5),1),'double');
    fwrite(fid,reshape(Sjk5,numel(Sjk5),1),'double');
    fwrite(fid,reshape(Skk5,numel(Skk5),1),'double');
    fclose(fid);

    fid = fopen(sprintf('%sS6.bin',DataOutput),'wb');
    fwrite(fid,[size(Sii6,1),size(Sii6,2),size(Sii6,3)],'uint16');
    fwrite(fid,reshape(Sii6,numel(Sii6),1),'double');
    fwrite(fid,reshape(Sij6,numel(Sij6),1),'double');
    fwrite(fid,reshape(Sik6,numel(Sik6),1),'double');
    fwrite(fid,reshape(Sjj6,numel(Sjj6),1),'double');
    fwrite(fid,reshape(Sjk6,numel(Sjk6),1),'double');
    fwrite(fid,reshape(Skk6,numel(Skk6),1),'double');
    fclose(fid);

    % fid = fopen(sprintf('%sS7.bin',DataOutput),'wb');
    % fwrite(fid,[size(Sii7,1),size(Sii7,2),size(Sii7,3)],'uint16');
    % fwrite(fid,reshape(Sii7,numel(Sii7),1),'double');
    % fwrite(fid,reshape(Sij7,numel(Sij7),1),'double');
    % fwrite(fid,reshape(Sik7,numel(Sik7),1),'double');
    % fwrite(fid,reshape(Sjj7,numel(Sjj7),1),'double');
    % fwrite(fid,reshape(Sjk7,numel(Sjk7),1),'double');
    % fwrite(fid,reshape(Skk7,numel(Skk7),1),'double');
    % fclose(fid);


    fid = fopen(sprintf('%sIJK2.bin',DataOutput),'wb');
    fwrite(fid,[size(SI2,1),size(SI2,2),size(SI2,3)],'uint16');
    fwrite(fid,reshape(SI2,numel(SI2),1),'uint16');
    fwrite(fid,reshape(SJ2,numel(SJ2),1),'uint16');
    fwrite(fid,reshape(SK2,numel(SK2),1),'uint16');
    fclose(fid);

    fid = fopen(sprintf('%sIJK3.bin',DataOutput),'wb');
    fwrite(fid,[size(SI3,1),size(SI3,2),size(SI3,3)],'uint16');
    fwrite(fid,reshape(SI3,numel(SI3),1),'uint16');
    fwrite(fid,reshape(SJ3,numel(SJ3),1),'uint16');
    fwrite(fid,reshape(SK3,numel(SK3),1),'uint16');
    fclose(fid);

    fid = fopen(sprintf('%sIJK4.bin',DataOutput),'wb');
    fwrite(fid,[size(SI4,1),size(SI4,2),size(SI4,3)],'uint16');
    fwrite(fid,reshape(SI4,numel(SI4),1),'uint16');
    fwrite(fid,reshape(SJ4,numel(SJ4),1),'uint16');
    fwrite(fid,reshape(SK4,numel(SK4),1),'uint16');
    fclose(fid);

    fid = fopen(sprintf('%sIJK5.bin',DataOutput),'wb');
    fwrite(fid,[size(SI5,1),size(SI5,2),size(SI5,3)],'uint16');
    fwrite(fid,reshape(SI5,numel(SI5),1),'uint16');
    fwrite(fid,reshape(SJ5,numel(SJ5),1),'uint16');
    fwrite(fid,reshape(SK5,numel(SK5),1),'uint16');
    fclose(fid);

    fid = fopen(sprintf('%sIJK6.bin',DataOutput),'wb');
    fwrite(fid,[size(SI6,1),size(SI6,2),size(SI6,3)],'uint16');
    fwrite(fid,reshape(SI6,numel(SI6),1),'uint16');
    fwrite(fid,reshape(SJ6,numel(SJ6),1),'uint16');
    fwrite(fid,reshape(SK6,numel(SK6),1),'uint16');
    fclose(fid);

    % fid = fopen(sprintf('%sIJK7.bin',DataOutput),'wb');
    % fwrite(fid,[size(SI7,1),size(SI7,2),size(SI7,3)],'uint16');
    % fwrite(fid,reshape(SI7,numel(SI7),1),'uint16');
    % fwrite(fid,reshape(SJ7,numel(SJ7),1),'uint16');
    % fwrite(fid,reshape(SK7,numel(SK7),1),'uint16');
    % fclose(fid);

end

if streamlines
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % ComputeStreamlines.m
    %
    % This script loads and processes structure tensor components at a
    % specified resolution to determine streamlines as a visible expression of
    % fibre direction.
    %
    % Updated by: Mark Trew, June 2022.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Set parameters and paths
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Level = 5; % frequency resolution of ST/Hessian data to use

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Data and path locations
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    InputPath = src_dir + '/binary/';
    OutputPath = src_dir + '/binary/';
    MaskPath = src_dir + '/mask/';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Load data
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(' ... loading in data ... \n');

    fid = fopen(sprintf('%sS%1d.bin',InputPath,Level),'r');
    N = fread(fid,3,'uint16');
    fprintf(' ... Dimension: [%d,%d,%d] ... \n',N);
    d2Xs = double(fread(fid,prod(N),'double'));
    dXYs = double(fread(fid,prod(N),'double'));
    dXZs = double(fread(fid,prod(N),'double'));
    d2Ys = double(fread(fid,prod(N),'double'));
    dYZs = double(fread(fid,prod(N),'double'));
    d2Zs = double(fread(fid,prod(N),'double'));
    fclose(fid);

    fid = fopen(sprintf('%sIJK%1d.bin',InputPath,Level),'r');
    N = fread(fid,3,'uint16');
    fprintf(' ... Dimension: [%d,%d,%d] ... \n',N);
    I = fread(fid,prod(N),'uint16');
    J = fread(fid,prod(N),'uint16');
    K = fread(fid,prod(N),'uint16');
    fclose(fid);

    % Read in mask data
    mask_paths = getImagePaths(MaskPath, extension);
    I3D = loadImageStack(mask_paths);
    [Nj, Ni, Nk] = size(I3D);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Manipulate loaded data
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I3D = permute(I3D,[2,1,3]);
    I3D = reshape(I3D,Ni*Nj*Nk,1);
    GD = ((K-1)*Nj+J-1)*Ni+I;
    MaskGD = find(I3D(GD));

    % Only use data within masked tissue
    d2Xs = d2Xs(MaskGD);
    dXYs = dXYs(MaskGD);
    dXZs = dXZs(MaskGD);
    d2Ys = d2Ys(MaskGD);
    dYZs = dYZs(MaskGD);
    d2Zs = d2Zs(MaskGD);
    I = I(MaskGD);
    J = J(MaskGD);
    K = K(MaskGD);

    % Index step sizes
    DJ = J(1+N(1))-J(1); DK = 225; %DK = K(1+N(1)*N(2))-K(1);
    DI = 32; % fix to this value regardless of data level.

    fprintf('DI: %d, DJ: %d, DK: %d\n',DI,DJ,DK);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Set up sample grid (seed points) for streamlines. The density of the seed
    % points is arbitrary.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('... Set up sample grid: ');

    [SI,SJ,SK] = ndgrid(1:round(2*DI/0.5):Ni,1:round(2*DI/0.5):Nj,1:round(2*DI/0.5):Nk);
    fprintf('%dX%dX%d\n',size(SI));

    NS = size(SI);
    LIS = sub2ind([Ni,Nj,Nk],reshape(SI,[prod(NS),1]),reshape(SJ,[prod(NS),1]),reshape(SK,[prod(NS),1]));
    IdxS = find(I3D(LIS));

    % Set up interpolant
    Fd2Xs = scatteredInterpolant(I,J,K,d2Xs);
    FdXYs = scatteredInterpolant(I,J,K,dXYs);
    FdXZs = scatteredInterpolant(I,J,K,dXZs);
    Fd2Ys = scatteredInterpolant(I,J,K,d2Ys);
    FdYZs = scatteredInterpolant(I,J,K,dYZs);
    Fd2Zs = scatteredInterpolant(I,J,K,d2Zs);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Calculate streamline paths
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine paths
    Paths = cell(length(IdxS),1);
    DS = 5;

    FiberIndex = 1; % fiber - otherwise use 1 for smallest eigenvalue/fiber
    MaxTrackLength = 10000; % Fiber tracks
    %MaxTrackLength = 500; % sheet tracks
    parfor i=1:length(IdxS)
        if ~mod(i,10) fprintf('Path: %d\n',i); end
        Paths{i} = FiberTrack([SI(IdxS(i)),SJ(IdxS(i)),SK(IdxS(i))],DS,I,J,K,Fd2Xs,FdXYs,FdXZs,Fd2Ys,FdYZs,Fd2Zs,I3D,[Ni,Nj,Nk],FiberIndex,MaxTrackLength);

    end

    if FiberIndex == 2
        exfname = sprintf('%s/StreamlinesSheet',OutputPath);
    elseif FiberIndex == 3
        exfname = sprintf('%s/StreamlinesNormal',OutputPath);
    else
        exfname = sprintf('%s/Streamlines_L%1d_FB',OutputPath,Level);
    end
    groupname = 'Streamlines';
    ExportStreamlines(Paths,[],exfname,groupname,1,1);
end

end