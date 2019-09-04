% Takes Function E and set inf&NaN values to 0
% then returns function
function E=handle_NaNInf(E)
E(isnan(E))=0;
E(isinf(E))=0;
end