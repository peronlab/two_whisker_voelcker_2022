%
% SP Jan 2011
% 
% This is a human-optimal Colormap, care of Dan O'Connor.
%
% USAGE:
%
%  cm = colormap_human(n)
%
%  n: how many elements you want in the colromap (default 64)
%
function cm = colormap_human(n)
  if (nargin < 1) ; n = 64 ; end

  % setup
	sig = round(0.1*n);
	np25 = round(0.25*n);
	np5 = round(0.5*n);
	np75 = round(0.75*n);

	x = 1:n; 
	r = normpdf(x,n,sig)+normpdf(x,np75,sig); r = r/max(r);
	g = normpdf(x,np5,sig)+normpdf(x,np75,sig)+normpdf(x,np25,sig); g = g/max(g);
	b = normpdf(x,0,sig)+normpdf(x,np25,sig); b = b/max(b);
	cm = [r' g' b'];

	% to see it alive
	% figure; plot(x,r,'r-',x,g,'g-',x,b,'b-'); xlim([1 n])

