function [bei, ie12, ie21] = binary_epsilon_indicator(pf1, pf2)
    ie12 = epsilon_indicator(pf1,pf2);
    ie21 = epsilon_indicator(pf2,pf1);
    bei = ( abs(ie12 - 1) < 1e-6 ) && ( (ie21 - 1) > 1e-6 );
end

function [min_factor] = epsilon_indicator(pf1, pf2)
    min_factor = max(bsxfun(@rdivide, max(abs(pf1),[],1), max(abs(pf2),[],1)));
end
