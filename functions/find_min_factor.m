function min_factor = find_min_factor(min_number, max_number)
facs = zeros(1, max_number - min_number);
fac_max = facs;
for index = min_number:max_number
    facs(index - min_number + 1) = length(factor(index));
    fac_max(index - min_number + 1) = max(factor(index));
end

min_factor = min(min_number + find(fac_max == min(fac_max)) - 1);

end