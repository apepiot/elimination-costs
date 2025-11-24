function y = myUnique(x)
  y = unique(x);
  y=y{:,:};
  if any(isnan(y))
    y(isnan(y)) = []; % remove all nans
    y(end+1) = NaN; % add the unique one.
  end
end