function B=GaussBlur1d(A,dw,dim)


    nt=size(A,dim);
    if rem(nt,2)==0
        kt=-nt/2+0.5:nt/2-0.5;
    else
        kt=-(nt-1)/2:(nt-1)/2;
    end
    if dim==1
        Fil=exp(-kt'.^2/dw^2)*ones(1,size(A,2));
    end
    if dim==2
        Fil=ones(size(A,1),1)*exp(-kt.^2/dw^2);
    end
    if dim==3
        for i=1:nt
            Fil(:,:,i)=ones(size(A,1),size(A,2))*exp(-kt(i).^2/dw^2);
        end
    end
    tfA=fftshift(fft(A,[],dim));
    B=real(ifft(ifftshift(tfA.*Fil),[],dim));
