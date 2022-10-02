function imgOut = img2rgb(img, Rvec, Gvec, Bvec)
    imgOutR = img(:, :, 1);
    imgOutG = img(:, :, 2);
    imgOutB = img(:, :, 3);
    Rvec = [0, Rvec, 256];
    Gvec = [0, Gvec, 256];
    Bvec = [0, Bvec, 256];

    for i = 1:size(Rvec, 2) - 1
        at = find(imgOutR(:, :) >= Rvec(i) & imgOutR(:, :) < Rvec(i + 1));
        imgOutR(at) = Rvec(i);
    end

    for i = 1:size(Gvec, 2) - 1
        at = find(imgOutG(:, :) >= Gvec(i) & imgOutG(:, :) < Gvec(i + 1));
        imgOutG(at) = Gvec(i);
    end

    for i = 1:size(Bvec, 2) - 1
        at = find(imgOutB(:, :) >= Bvec(i) & imgOutB(:, :) < Bvec(i + 1));
        imgOutB(at) = Bvec(i);
    end

    imgOut(:, :, 1) = imgOutR;
    imgOut(:, :, 2) = imgOutG;
    imgOut(:, :, 3) = imgOutB;
end
