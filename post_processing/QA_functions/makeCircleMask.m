function [circleMask] = makeCircleMask(im_size, Xc, Yc, Radius)
circleMask = zeros(im_size) ;
I = im_size(1) ; J = im_size(2) ;
for i = 1:I
    for j = 1:J
        distance = (double(i) - Xc).^2 + (double(j) - Yc).^2 ;
        if distance <= Radius.^2
            circleMask(i, j) = 1 ;
        end
    end
end
end