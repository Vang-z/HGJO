function PSNRV = PSNR(Iin, Iout)
    Iin = double(Iin);
    Iout = double(Iout);
    
    [M, N] = size(Iin);
    error = Iin - Iout;
    MSE = sum(sum(error .* error)) / (M * N);
    if(MSE > 0)
        PSNRV = 10 * log(255^2 / MSE) / log(10);
    else
        PSNRV = 99;
    end
    PSNRV = PSNRV(:)';
end
