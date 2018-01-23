function normalizeRoiTrace2WindowedRange(R)

 X = [R.Trace];
   fs=20;
   winsize = 1*fs;
   numwin = floor(size(X,1)/winsize)-1;
   xRange = zeros(numwin,size(X,2));
   for k=1:numwin
	  windex = (winsize*(k-1)+1):(winsize*(k-1)+20);
	  xRange(k,:) = range(detrend(X(windex,:)), 1);
   end
   X = bsxfun(@rdivide, X, mean(xRange,1));
   for k=1:numel(R)
	  R(k).Trace = X(:,k);
   end


