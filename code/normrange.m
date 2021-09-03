function X = normrange(X)
%NORMRANGE Normalize data by mapping its range to [0,1]


xmax = max(X);
xmin = min(X);

X=(X-xmin)/(xmax-xmin);
end

