%�G�[�W�F���g��ʂ��Ƃ�F���L�q
%x:[n*1]agt_type:[n*1]factor:[n*num_factor]
function outF = eF(x,A,B)

outF = 0;
n = numel(x);
for i = 1:n
	outF=outF+A(i)*(x(i)-B(i))^2;
end