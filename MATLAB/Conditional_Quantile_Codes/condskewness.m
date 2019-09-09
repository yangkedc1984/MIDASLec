function [cskew] = condskewness(condup,conddown,condmed,qlevel)
%Conditional Skewness function
  num = condup + conddown - 2 * condmed;
  den = condup - conddown;
  %Kornish-Fisher constant:
  const = 6/norminv(qlevel);
  cskew = num./den*const;
end

