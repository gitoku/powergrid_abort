%エージェント種別ごとにdF/dxiを記述
%x:[1]agt_type:[1]factor:[num_factor]
function out = dFdx(agt_type,x,factor)

	cA = factor(1);
    ca = factor(2);
    cb = factor(3);
    
	sa = factor(4);
	sb = factor(5); %影響なし
    sc = factor(6);
    
	da = factor(7);
	db = factor(8); %影響なし
	dc = factor(9);

%エージェント種別###
%1:需要家
%2:供給家
%3:送電家


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

%下に凸に変換
out = -out;

% if out > .9
%    disp('er')
%     pause on
% else 
%     
% disp([out agt_type a b c x]);
% end