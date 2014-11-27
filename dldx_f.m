function out = dFdx(agt_type,x,factor)

if agt_type_i == 1
	out = factor(1) * ( factor(2) * x + factor(3) )^2
elseif agt_type_i == 2
	out = factor(1) * ( factor(2) * x - factor(3) )^2
elseif agt_type_i == 3
	out = factor(1) * ( factor(2) * x + factor(3) - factor(4) )^2
end