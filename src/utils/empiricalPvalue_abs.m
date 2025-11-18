function pvalue = empiricalPvalue_abs(empiricalDist, observedValue)

arguments
    empiricalDist (:,:) double
    observedValue (1,1) double
end

N = numel(empiricalDist);
empiricalDist_abs = abs(empiricalDist);
observedValue_abs = abs(observedValue);
k = sum(empiricalDist_abs >= observedValue_abs);
pvalue = (k+1) / (N+1);

end
