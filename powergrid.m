
%% 初期化処理
f_init = input('initialize?[y,n]','s');   %'y'の時,初期化処理を行う
if isempty(f_init)
    f_init = 'n';
end


if f_init=='y' || exist('define.mat','file')~=2   %'define.mat'がなければ強制的に初期化処理
    clear all;
    disp('初期化中');

    %% グラフの定義###
    
    % エッジ集合にて定義
   E = [
        1,4;2,5;3,6;4,5;4,7;5,6;6,9;6,34;7,8;8,9;8,11;9,12;9,34;10,11;11,12;12,39;...   %region1
        13,14;14,15;14,17;14,34;15,16;16,19;16,35;17,18;17,20;17,18;18,19;19,37;20,38;... %region2
        21,22;21,35;22,23;22,24;23,25;24,25;24,36;25,36;... %region3
        26,36;26,37;... %region4
        27,38;27,39;27,40;... %region5
        28,29;28,31;29,30;29,32;30,33;30,40;32,33;33,40; %region6
    ];

%   小規模グラフ
%   E = {{1,4},{2,3},{3,4},{3,5},{3,7},{3,9},{4,9},{4,15},{6,8},{7,8},{7,10},{7,12},{8,9},{8,12},{9,14},{11,12},{12,13},{15,16},{15,17}};

    num_x = max( max(E) ); %xの要素数


    % 隣接行列 N
    for i=1:length(E);
        N( E(i,1) , E(i,2) ) = 1;
        N( E(i,2) , E(i,1) ) = 1;  %%無向グラフの場合
    end


    %% グラフラプラシアン L
    L_diag=zeros([1 num_x]);
    for i=1:num_x
        for j=1:length(E)
%             if E(j,1)==i   %有向グラフ
            if E(j,1)==i || E(j,2)==i   %無向グラフ
                L_diag(i) = L_diag(i)+1;
            end
        end
    end

    Lp = -N;
    for i=1:num_x
        Lp(i,i)=L_diag(i);
    end



    %% xの決定
    x = sym('x',[num_x 1]);

    %エージェント種別###
    %1:供給家[3:solar,4:wind]
    %2:需要家[1:home,2:factory]
    %3:送電家
   agt_type=[
        1,1,1,2,2,2,2,2,2,2,2,2,...%region1
        1,2,2,2,2,2,2,2,...%region2
        2,2,2,2,2,...%region3
        1,...%region4
        1,...%region5
        2,2,2,2,2,2,...%region6
        3,3,3,3,3,3,3%送電家
    ];
    agt_sub_type = [
        3,5,5,1,1,1,1,1,1,1,1,1,...
        3,1,1,1,1,1,1,1,...
        2,2,2,2,2,...
        3,...
        4,...
        1,1,1,1,1,1
    ];

% 小規模グラフ
% agt_type =      [1 2 3 3 1 2 3 3 3 1 2 3 1 1 3 1 2];
% agt_sub_type =  [3 1 0 0 3 1 0 0 0 3 2 0 3 4 0 3 1];



    %%  Gの決定###
    G_hat_sym = sym('G_hat_sym',[num_x 1]);
    G_hat_sym = [
        x(1)+x(2)+x(3)-x(4)-x(5)-x(6)-x(7)-x(8)-x(9)-x(10)-x(11)-x(12)+x(34)+x(39);
        x(13)-x(14)-x(15)-x(16)-x(17)-x(18)-x(19)-x(20)-x(34)+x(35)+x(37)+x(38);
        -x(21)-x(22)-x(23)-x(24)-x(25)-x(35)+x(36);
        x(26)-x(36)-x(37);
        x(27)-x(38)-x(39)+x(40);
        -x(28)-x(29)-x(30)-x(31)-x(32)-x(33)-x(40);
    ];
    num_G = length(G_hat_sym)*2;

%       小規模グラフ
%     G_hat_sym = [
%         x(1)-x(2)-x(3)+x(4);
%         x(3)+x(5)-x(6)+x(7)+x(8)+x(9);
%         -x(7)+x(10)-x(11)+x(12);
%         -x(8)-x(12)+x(13);
%         -x(4)-x(9)+x(14)-x(15);
%         x(15)+x(16)-x(17)
%         ];
    
    
    % 初回
    if exist('G.mat','file')~=2
        disp('Gの準備中');
        G_sym = [G_hat_sym.' (-G_hat_sym).'].';

        %Gをxの要素ごとに分解(Gmatrixの列に対応)
        Gmatrix_sym = sym('Gmatrix',[num_G num_x]);
        for i=1:num_x
            Gmatrix_sym(:,i) = G_sym;
        end
        for l=1:num_G
            for i=1:num_x
               for j=1:num_x
                   if i~=j
                        Gmatrix_sym(l,i) = subs(Gmatrix_sym(l,i),x(j),0);
                   end
               end
            end
        end

       for m=1:num_G
            for n=1:num_x
                Gm{m,n} = matlabFunction(Gmatrix_sym(m,n),'vars',{x(n)});
            end
            G{m} = matlabFunction(G_sym(m),'vars',{x});
       end
       save('G','G_sym','Gmatrix_sym','Gm','G');

    % 二回目以降
    else
        load('G');
    end
    G_sym
    
    
    %% λの定義
    num_lambda = num_G;
    lambda = sym('lambda',[num_lambda 1]);







    %% λGの決定
    lG_sym = lambda.'*G_sym;

    %%  d/dx (λG)の決定
    dlGdx_sym = sym('dlGdx_sym',[num_x 1]);
    for n=1:numel(x)
       dlGdx_sym(n) = diff(lG_sym,x(n));
       dlGdxi{n} = matlabFunction(dlGdx_sym(n),'vars',{lambda});
    end
    dlGdx = matlabFunction(dlGdx_sym,'vars',{lambda});


    
    %%  Region定義
    lambda_matrix = zeros([num_x num_lambda]);
    for i=1:num_x
        for j=1:num_lambda
            if isnan(subs(G_sym(j),x(i),NaN))
                lambda_matrix(i,j) = 1;
            end
        end
    end

    % Niの導出 (現状使用していない)
    v=lambda_matrix.';
    Ni=cell([num_lambda/2 1]);
    for i=1:num_lambda/2
        size = sum(v(i,:));
        Ni(i)={zeros(size)};
        cntj=1;
        cntk=1;
        for j=1:num_x
            for k=1:num_x
                if v(i,j)==1 && v(i,k)==1
                    Ni{i}(cntj,cntk)=N(j,k);
                    %disp([cntj cntk Lp(j,k)]);
                    cntk=cntk+1;
                    if cntk>size
                        cntj=cntj+1;
                        cntk=1;
                    end
                end
            end
        end
    end

    save('define','num_x','num_lambda','N','Lp','L_diag',...
        'lambda_matrix','agt_type','agt_sub_type',...
        'G_sym','Gmatrix_sym','Gm','G','dlGdxi','dlGdx_sym','lG_sym'...
        );
    clear all;
    disp('初期化完了');
end
clear f_init;

clear all









load('define');
rng(3564);


%% パラメータ設定###
A =10;
B = .02;
gamma = .0005;
c = .1./L_diag;

B_p = B/sum(1./c);  %スーパバイザ用のB


day = 24;
c_delay = 0;
stp_max = day*4+1;    %s(実行step数)の最大
eps_x = .001;   %x[k]の更新の打ち切り基準:dx[k]<eps_x
eps_t = .001;    %θ[k]の更新の打ち切り基準:[{max(θ[k])-min(θ[k])}/mean(θ[k])]<eps_t
dx_max = 1000;    %x[k]の更新の計算中止dx
kt_max = 1000;     %θ[k]の更新の計算打ち切りk

%θの合意の経過データ用###
wtc_m = 2;      %監視するθiのi
wtc_step = 2;   %監視するθiのs
g = rand([num_x,1])*1;


%% シミュレーション実行
f_run = input('run?[y,n]','s');  %'y'で実行
if isempty(f_run)
    f_run = 'n';
end
if f_run == 'y'
    %% x,λの推移を記憶
    % 
    % $$e^{\pi i} + 1 = 0$$
    % 
    X = ones(num_x,stp_max);
    X_min = ones(num_x,stp_max*60);
    LAMBDA = cell(num_x,stp_max);
    LAMBDA_s = cell(num_lambda,stp_max);    % λ検算用(スーパバイザ方式)

    %% 初期条件(step = 1)
    %rand('seed',100);%乱数固定用
    X(:,1) = rand([num_x,1]);   %[0~1]の乱数
    for i=1:num_x
        for mi=1:60
            X_min(:,mi) = X(:,1);
        end
        
        for j=1:num_lambda
          LAMBDA{i,1}(j,1) = lambda_matrix(i,j)*rand(1); %[0~1]の乱数
        end
    end

    for j=1:num_lambda
        SUM = 0;
        for i=1:num_x
            SUM = SUM + LAMBDA{i,1}(j); % λ検算用(スーパバイザ方式)
        end
        LAMBDA_s(j,1) = {SUM/sum(lambda_matrix(:,j))};
    end

    

    
    %% ステップ実行(step >= 2)
    disp('実行中...')
    for step = 2:stp_max
    
        

            
        %% xの更新
        %x[0]の準備<(17)_a>
        x = X(:,step-1);
        for i=1:num_x   %各ノードについて
            % x[k]がε以下となるまで更新<(17)_b>
            kx=0;
            while kx < 60
                if step==2
                    l = LAMBDA{i,step-1};
                elseif kx<=c_delay
                    l = LAMBDA{i,step-2};
                else
                    l = LAMBDA{i,step-1};
                end
                    

                f = 1;
                
                now = (step-1)*60+kx;

                if agt_type(i)==2
                    if agt_sub_type(i)==1
%                         f = home1(now,day*60,2,1,3);
                        f = 2;
                    elseif agt_sub_type(i)==2
                        f = 2.5;
                    end
                elseif agt_type(i)==1
                    if agt_sub_type(i)==3
%                         f = solar(now,day*60,5);
                        f=3;
                    elseif agt_sub_type(i)==4
%                         f = wind(now,day*60,3,2);
                        f=2;
                    elseif agt_sub_type(i)==5
%                         f = home1(now,day*60,17,6,10);   %6,1,6
                        f=10;
                    end
                end


                factor = [g(i) f 0  f 0 g(i)*2  5 0 2];
                df = dFdx( agt_type(i),  x(i), factor(:) );
                dg = dlGdxi{i}(l);
                x(i) = x(i) - A* ( gamma*df + dg );
                

                kx=kx+1;
                
                X_min(i,(step-1)*60+kx) = x(i);
            end
        end
        %xの更新<(17)_c>
        X(:,step) = x;





         %% λの更新
        % λ検算用(スーパバイザ方式)
        for m = 1:num_lambda
            LAMBDA_s(m,step) = {max(0, LAMBDA_s{m,step-1}+B_p*G{m}( X(:,step)))};
        end

        theta = zeros([num_x num_lambda]);
        next_theta = zeros([num_x num_lambda]);

        % θi[0]の準備<(28)_a>
        for n=1:num_x
            for m=1:num_lambda
                if lambda_matrix(n,m)
                    theta(n,m) = LAMBDA{n,step-1}(m) + c(n)*B*Gm{m,n}(X(n,step));
                else
                    theta(n,m)=NaN;
                end
            end
        end
        
        
        
        % θi[k]の相対差がε以下となるまで更新<(28)_b>
        for m=1:num_lambda
            kt=1;
            while true
                for n = 1:num_x
                    sigma=0;
                    for o = 1:num_x
                        if (N(n,o)==1 && ~isnan(theta(o,m)) && n~=o)
                            sigma = sigma + lambda_matrix(n,m)*(theta(n,m)-theta(o,m));
                        end
                    end
                    next_theta(n,m) = theta(n,m) - c(n)*sigma;
                end
                theta(:,m) = next_theta(:,m);

                converge = abs((max(theta(:,m))-min(theta(:,m)))/nanmean(theta(:,m)));

                if converge < eps_t     %合意達成
                    break;
                elseif kt>kt_max
                    disp('θ計算打ち切り');
                    disp([step m converge]);
                    break;
                end

                % 特定のsについてのθ(m番目要素)を記録
                if m==wtc_m && step==wtc_step
                    THETA(:,kt) = theta(:,m);
                    KT_END = kt;
                end
                kt=kt+1;
            end
        end

        %λの推定値の更新<(28)_c>
        for n = 1:num_x
            LAMBDA(n,step) = {theta(n,:).'};
        end
        

        % 経過の表示
        % s,k(xについて),k(thetaについて)
        disp([step kx kt]);
    end

    save('result','stp_max','X','LAMBDA','LAMBDA_s','KT_END','THETA','wtc_m','wtc_step','c','factor','X_min','day');
    clear all;
end
clear f_run;









%% 結果の表示
f_plot = input('plot?[y,n]','s');  %'y'で実行
if isempty(f_plot)
    f_plot = 'n';
end
if f_plot == 'y'
    load('define');
    load('result');

    % θの合意の経過についてのデータ整形
    k_max = floor((KT_END-1)/10)*10;  %グラフ右端が10の倍数になるようにデータ切り捨て
    dec = k_max+100;
%     k = 0:k_max;
k=0:dec-1;
    theta = zeros([num_x k_max+1]);
    theta_s = zeros([k_max+1]);
    for kt = 1:k_max+1
        theta(:,kt) = THETA(:,kt);
        theta_s(kt) = THETA(1,kt); %スーパバイザ方式の結果(インチキ)
    end
    
    for kt = k_max+2:dec
        theta(:,kt) = theta(:,k_max+1);
        theta_s(kt) = theta(1,k_max+1); %スーパバイザ方式の結果(インチキ)
    end


    LAMBDA_plot = zeros(num_lambda,stp_max);
    for t=1:stp_max
        for i=1:num_lambda
            LAMBDA_plot(i,t) = LAMBDA_s{i,t};
        end
    end
    
    LAMBDA_min = zeros([num_lambda stp_max*60]);
    for t=1:stp_max
        for m=1:60
            LAMBDA_min(:,(t-1)*60+m) = LAMBDA_plot(:,t);
        end
    end

    % F,Gの推移についてのデータ計算
    FX = zeros([stp_max 1]);
    GX = zeros([num_lambda stp_max]);
    for step = 1:stp_max
        % FX(step) = F(X(:,step));
        FX(step) = eF(X(:,step),agt_type,factor);
        for m=1:num_lambda
            GX(m,step) = G{m}(X(:,step));
        end
    end
    
    
    GX_min = zeros([num_lambda stp_max*60]);
    for m=1:num_lambda
        for step=1:stp_max*60
            GX_min(m,step) = G{m}(X_min(:,step));
        end
    end
    
%     LAMBDA_min = zeros([num_lambda stp_max]);
%     for step = 1:stp_max
%         for m=1:num_lambda
%             for step=1:stp_max*60
%                 GX_min(m,step) = G{m}(X_min(:,step));
%             end
%         end
%     end

    ofs=1;
    time = -ofs:stp_max-1-ofs;  %s
    
    mmax = stp_max*60;
    time_min = 1:mmax;
    
    time_h = time_min./60;

%     figure(1);
%     plot(time,X,'LineWidth',1.5);
%     title('x');
%     grid on;
%     xlim([0 stp_max-1-ofs]);
%     set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;

%     figure(2);
%     plot(time,LAMBDA_plot,'LineWidth',1.5);
%     title('λ');
%     grid on;
%     xlim([0 stp_max-1-ofs]);
%     set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;

%     figure(3);
%     plot(time,FX,'LineWidth',1.5);
%     %title('F');
%     grid on;
%     xlim([0 stp_max-1-ofs]);
%     set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
% 
%     figure(4);
%     plot(time,GX,'LineWidth',1.5);
%     title('G');
%     grid on;
%     xlim([0 stp_max-1-ofs]);
%     set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;


%     legend('1','2','3','4','5','6');
    
%         figure(11);
%     plot(time_min,LAMBDA_min,'LineWidth',1.5);
%     title('λ');
%     grid on;
%     set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
%     legend('1','2','3','4','5','6');
    
    
%     figure(6);
%      plot(time_min,X_min(1:6,:),'LineWidth',1.5);
%     title('x 1');
%     grid on;
%     set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
%     legend('1+s','2-h','3=','4=','5+s','6-h');
%     
%      figure(8);
%      plot(time_min,X_min(7:11,:),'LineWidth',1.5);
%     title('x 2');
%     grid on;
%     set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
%     legend('7=','8=','9=','10-f','11+f');
%     
%     figure(9);
%      plot(time_min,X_min(12:17,:),'LineWidth',1.5);
%     title('x 3');
%     grid on;
%     set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
%     legend('12=','13+s','14+w','15=','16+s','17-h');
%     

    
%     figure(12);
%     plot(time_min,GX_min(7:12,:),'LineWidth',1.5);
%     title('G-');
%     set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
%     grid on;
%     legend('1','2','3','4','5','6');



    figure(20);
    plot(time_h,X_min(list,:),'LineWidth',1.5);
    title('region1');
    grid on;
    set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
    legend('1+','2+','3+','4-','5-','6-','7-','8-','9-','10-','11-','12-','34=+','39=+');
    
    
 
    
    
    figure(1);
    list = [1 2 3 4 5 6 7 8 9 10 11 12 34 39];
    plot(time_h,X_min(list,:),'LineWidth',1.5);
    title('region1');
    grid on;
    set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
    legend('1+','2+','3+','4-','5-','6-','7-','8-','9-','10-','11-','12-','34=+','39=+');
    
    
        figure(2);
        list = [13 14 15 16 17 18 19 20 34 35 37 38];
     plot(time_h,X_min(list,:),'LineWidth',1.5);
    title('region2');
    grid on;
    set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
    legend('13+','14-','15-','16-','17-','18-','19-','20-','34=-','35=+','37=+','38=+');
    
    figure(3);
    list = [21 22 23 24 25 35 36];
    plot(time_h,X_min(list,:),'LineWidth',1.5);
    title('region3');
    grid on;
    set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
    legend('21-','22-','23-','24-','25-','35=-','36=+');
    
        figure(4);
    list = [26 36 37];
    plot(time_h,X_min(list,:),'LineWidth',1.5);
    title('region4');
    grid on;
    set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
    legend('26+','36=-','37=-');
    
        figure(5);
    list = [27 38 39 40];
    plot(time_h,X_min(list,:),'LineWidth',1.5);
    title('region5');
    grid on;
    set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
    legend('27+','38=-','39=-','40=+');
    
        figure(6);
    list = [28 29 30 31 32 33 40];
    plot(time_h,X_min(list,:),'LineWidth',1.5);
    title('region6');
    grid on;
    set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
    legend('28-','29-','30-','31-','32-','33-','40=-');

    save('result2','GX','FX','LAMBDA_plot');
    
    LAMBDA2 = LAMBDA_min(1:6,:)-LAMBDA_min(7:12,:);
    
    
        figure(8);
    plot(time_h,GX_min(1:6,:),'LineWidth',1.5);
%     title('G+(余ると上昇)');
    set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
    xlim([0 70]);
    grid on;
%     legend('1','2','3','4','5','6');
    
    figure(7);
    plot(time_h,LAMBDA_min(7:12,:)-LAMBDA_min(1:6,:),'LineWidth',1.5);
%     title('λ- - λ+(不足で上昇)');
    grid on;
    set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
    xlim([0 70]);
    
       figure(10);
    kh = k.*(0.5/dec);
    plot(kh,theta,'LineWidth',1.5);
    axis 'auto y';
    hold on;
%     plot(k,theta_s,'--','LineWidth',1.5);
    hold off;
    grid on;
    axis([0 0.5 0.1 .7]);
    set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
    
%     
%     f =zeros([1 (stp_max)*60]);
%     tt = 1:(stp_max)*60;
%     for i=1:(stp_max)*60
%         f(i) = home1(tt(i),day*60,2,1,3);
%     end
%     
%     tth = tt./60;
%     figure(11);
%     plot(tth,f,tth,X_min(5:12,:),'LineWidth',1.5);
%     grid on;
%     xlim([0 70]);
%     legend('h','5','6','7','8','9','10','11','12');
%     set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;
    
    
%         figure(12);
%     plot(tth,f,tth,X_min(8,:),'LineWidth',1.5);
%     grid on;
%     xlim([0 70]);
%     set(gca,'FontName','Times','FontSize',18,'LineWidth',1.5) ;

end

% f_export = input('export?[number]','s');
% if isempty(f_export)
%     f_export = 'n';
% end
% if f_export ~= 'n'
%     out=[X_min.' GX_min.' LAMBDA2.'];
%     name=strcat('data',f_export,'.csv');
%     csvwrite(name,out);
% end

clear f_plot;