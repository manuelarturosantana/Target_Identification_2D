function [dists, inds] = vdist(v1,v2)
    % Function which computes the index of the closest element of v1 to v2
    v1 = v1(:); v2 = v2(:);
    dists = []; inds = [];
    for ii = 1:length(v1)
        [di,ind] = min(abs(v1(ii) - v2));
        dists(ii) = di;
        inds(ii) = ind;
    end
end