
function [st] = genpsfstruct(f,pixelsize,dz,datatype)
% calculate f,dfx,dfy,dfxy,dfz,dfxz,dfyz,dfxyz at each sample point
N = size(f,1);
Nz = size(f,3);
dv = pixelsize;
dfx = zeros(N,N,Nz);
dfy = zeros(N,N,Nz);

dfx(:,1,:) = (f(:,2,:)-f(:,1,:))./dv;
dfx(:,N,:) = (f(:,N,:)-f(:,N-1,:))./dv;
dfy(1,:,:) = (f(2,:,:)-f(1,:,:))./dv;
dfy(N,:,:) = (f(:,N,:)-f(:,N-1,:))./dv;

for ii = 2:N-1
    dfx(:,ii,:) = (f(:,ii+1,:)-f(:,ii-1,:))./2./dv;
    dfy(ii,:,:) = (f(ii+1,:,:)-f(ii-1,:,:))./2./dv;
end


dfxy = zeros(N,N,Nz);
dfxy(1,:,:) = (dfx(2,:,:)-dfx(1,:,:))./dv;
dfxy(N,:,:) = (dfx(N,:,:)-dfx(N-1,:,:))./dv;

for ii = 2:N-1
    dfxy(ii,:,:) = (dfx(ii+1,:,:)-dfx(ii-1,:,:))./2./dv;
end


dfz = zeros(N,N,Nz);
dfxz = zeros(N,N,Nz);
dfyz = zeros(N,N,Nz);
dfxyz = zeros(N,N,Nz);

dfz(:,:,1) = (f(:,:,2)-f(:,:,1))./dz;
dfz(:,:,Nz) = (f(:,:,Nz)-f(:,:,Nz-1))./dz;
dfxz(:,:,1) = (dfx(:,:,2)-dfx(:,:,1))./dz;
dfxz(:,:,Nz) = (dfx(:,:,Nz)-dfx(:,:,Nz-1))./dz;
dfyz(:,:,1) = (dfy(:,:,2)-dfy(:,:,1))./dz;
dfyz(:,:,Nz) = (dfy(:,:,Nz)-dfy(:,:,Nz-1))./dz;
dfxyz(:,:,1) = (dfxy(:,:,2)-dfxy(:,:,1))./dz;
dfxyz(:,:,Nz) = (dfxy(:,:,Nz)-dfxy(:,:,Nz-1))./dz;


for ii = 2:Nz-1
    dfz(:,:,ii) = (f(:,:,ii+1)-f(:,:,ii-1))./2./dz;
    dfxz(:,:,ii) = (dfx(:,:,ii+1)-dfx(:,:,ii-1))./2./dz;
    dfyz(:,:,ii) = (dfy(:,:,ii+1)-dfy(:,:,ii-1))./2./dz;
    dfxyz(:,:,ii) = (dfxy(:,:,ii+1)-dfxy(:,:,ii-1))./2./dz;
end

switch datatype
    case 'vector'
        st.F = f(:);
        st.Fx = dfx(:);
        st.Fy = dfy(:);
        st.Fxy = dfxy(:);
        st.Fz = dfz(:);
        st.Fxz = dfxz(:);
        st.Fyz = dfyz(:);
        st.Fxyz = dfxyz(:);
    case 'matrix'
        st.F = f;
        st.Fx = dfx;
        st.Fy = dfy;
        st.Fz = dfz;
        st.Fxy = dfxy;
        st.Fxz = dfxz;
        st.Fyz = dfyz;
        st.Fxyz = dfxyz;
end



