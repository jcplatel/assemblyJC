function B=GaussBlur(A,dw)
%2D Gaussian filter in Fourier domain on 2 first dimensions
[nx,ny,nz]=size(A);
if rem(nx,2)==0
    kx=-nx/2+0.5:nx/2-0.5;
else
    kx=-(nx-1)/2:(nx-1)/2;
end
if rem(ny,2)==0
    ky=-ny/2+0.5:ny/2-0.5;
else
    ky=-(ny-1)/2:(ny-1)/2;
end
[Ky,Kx]=meshgrid(ky,kx);
Krho=sqrt(Kx.^2+Ky.^2);
Fil=exp(-Krho.^2/dw.^2);
B=zeros(nx,ny,nz);
tfA=zeros(nx,ny,nz);
parfor i=1:nz
    tfA(:,:,i)=fftshift(fft2(A(:,:,i)));
%     B(:,:,i)=abs(ifft2(ifftshift(tfA(:,:,i).*Fil)));
    B(:,:,i)=real(ifft2(ifftshift(tfA(:,:,i).*Fil)));
end
