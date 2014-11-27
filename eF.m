%�G�[�W�F���g��ʂ��Ƃ�F���L�q
%x:[n*1]agt_type:[n*1]factor:[n*num_factor]
function outF = eF(x,agt_type,factor)

	a = factor(1);
	b = factor(2);
	c = factor(3);

%�G�[�W�F���g���###
%1:���v��
%2:������
%3:���d��
outF = 0;
n = numel(x);
for i = 1:n
	if agt_type(i) == 1
		A = -(b/a - c) / a;
		B = a - b/(A*a);
% 		if x<a
% 			out = A*x*(x - B);
% 		else
% 			out = c*(x - a) + b;
%         end
        out = factor(1)*(factor(2)*x(i)-factor(3))^2;
	elseif agt_type(i) == 2
		out = factor(1)*(factor(2)*x(i)-factor(3))^2;
	elseif agt_type(i) == 3
		out = factor(1)*(factor(2)*x(i)+factor(3)-factor(4))^2;
	end
	outF=outF+out;
end