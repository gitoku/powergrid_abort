%�G�[�W�F���g��ʂ��Ƃ�dF/dxi���L�q
%x:[1]agt_type:[1]factor:[num_factor]
function out = dFdx(agt_type,x,factor)

	cA = factor(1);
    ca = factor(2);
    cb = factor(3);
    
	sa = factor(4);
	sb = factor(5); %�e���Ȃ�
    sc = factor(6);
    
	da = factor(7);
	db = factor(8); %�e���Ȃ�
	dc = factor(9);

%�G�[�W�F���g���###
%1:���v��
%2:������
%3:���d��


if agt_type == 2
 	if x<ca
        out=cA*(-2*x+2*ca);
    else
        out=cA*(.5*(-2*x+2*ca));
    end
  
elseif agt_type == 1
 	if x<sa
        out=sc*(-2*x+2*sa);
    else
        out=sc*10*(-2*x+2*sa);
    end
    
elseif agt_type == 3
 	if x<(-da)
        out=-dc*(2*x+2*da);
    elseif x<da
        out=0;
    else
        out=-dc*(2*x-2*da);
    end
end

%���ɓʂɕϊ�
out = -out;

% if out > .9
%    disp('er')
%     pause on
% else 
%     
% disp([out agt_type a b c x]);
% end