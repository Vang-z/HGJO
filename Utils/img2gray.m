function imgOut = img2gray(img, Rvec)
    limites = [0, Rvec, 255];
    tamanho = size(img);
    imgOut(:, :) = img * 0;
    for i = 1:tamanho(1, 1)
        for j = 1:tamanho(1, 2)
            for k = 1:size(limites, 2) - 1
                if(img(i, j) >= limites(1, k) && img(i, j) <= limites(1, k + 1))
                    imgOut(i, j, 1) = limites(1, k);
                end
            end
        end
    end
end
