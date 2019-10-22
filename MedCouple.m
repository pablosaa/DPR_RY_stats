function [MC, thLow, thHig] = MedCouple(data, varargin),

	% Function to calculate the MedCouple from a given data variable:
	% Based on pseudo-code from article in Wikipedia.
	%
	% (c) 2019, Pablo Saavedra G.
	% Geophysical Institue, Univesitu of Bergen

	if isempty(data),
		MC = NaN;
		thLow = -Inf;
		thHig = +Inf;
		return;
	end
	

	% Sorting in decreasing order can be done in-place in O(n log n) time
	X = sort(data, 'descend');
	xm = nanmedian(X);
	xscale = 2*max(abs(X));

	if nargin == 2;
		disp('factor assigned!')
		Wfaktor = varargin{1};
	else
		Wfaktor = 1.5;
	end
	IQR = quantile(data, [.25 .75]);
	
	% define the upper and lower centred and rescaled vectors
  % they inherit X's own decreasing sorting
	xminus = X(X<xm);
	xplus = X(X>xm);
	Zminus = (xminus-xm)/xscale;
	Zplus = (xplus-xm)/xscale;

	p=length(Zplus);
	q=length(Zminus);


	% define the kernel function closing over Zplus and Zminus
	H_ij = arrayfun(@(a,b) (a+b)/(a-b), Zplus, Zminus);

	% signum function not implemented for (i,j) | xplus_i == xm == xminus_j 
	%sign(p-1-i-j);
	
	MC = median(H_ij);

	if MC >= 0,
		thLow = IQR(1) - Wfaktor*exp(-4*MC)*diff(IQR);
		thHig = IQR(2) + Wfaktor*exp(-3*MC)*diff(IQR);
	elseif MC<0
		thLow = IQR(1) - Wfaktor*exp(-3*MC)*diff(IQR);
		thHig = IQR(2) + Wfaktor*exp(-4*MC)*diff(IQR);
	else
		error('MC non defined numeric value!');
	end
	
end
% end of function
