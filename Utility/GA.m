% 完成していません。

function [x,J0,LOG] = GA(func,x_dim,option)
% Solving  x_opt = argmin_x func(x) 
% xmin  <-  [xmin(1,1), ... , xmin(1,x_dim)]
% xは行ベクトル
% bitmask = true -> random
    arguments
        func
        x_dim {mustBeInteger}
        option.xmax     =  1e+4
        option.xmin     = -1e+4
        option.format   = "uint32"
        option.bitmask  = [0x0000ffff,0x00ffff00,0x0f0f0f0f]
    end

    % Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xbounds = [option.xmin; option.xmax];               % 想定領域
    Scale   = intmax(option.format);                    % 決定変数の階調

    N       = 100;                                      % 個体数
    G       = 10;                                       % 世代数
    GroupA  = (xbounds(2,:)-xbounds(1,:)).*rand(N,x_dim) + xbounds(1,:);
                                                        % 現世代
    Gain    = 1;                                        % 適応度計算パラメータ
    Bias    = 0;                                        % 適応度計算パラメータ
    
    % 交叉設定
    Parents_num = 20;                                   % ペアの数
    parents_selection = @(t)Roulette_wheel...
                (t,Gain,Bias,N,2*N,Parents_num);        % 親の選択関数
    bitmask_num = length(option.bitmask);               % default = 3　マスクの個数
    if isscalar(bitmask_num)
        bitmask = @(x)randi(Scale+1)-1;                 % ランダムなマスク作成
    else
        bitmask = option.bitmask;
    end

    N_cross = Parents_num*bitmask_num;                  % 交叉個数
    N_mutat = 16;                                       % 変異個数
    EliteCount = 3;                                     % エリートの個数

    % 元の集団A は N                                 個 100
    % 交叉集団Q1は N_cross(=bitmask_num*Parents_num) 個 60
    % 突然変異Q2は N_mutation                        個 16
    % 選択では N 個の選択を行う．
    % 全体の個体数を可変にする場合はRoulette_wheel等の引数を増やす必要あり
    
    % Initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Calc_Encode = @(x)Encode(x,xbounds,option.format,Scale);
    Calc_Decode = @(x)Decode(x,xbounds,Scale);
    Crossover   = @(x,s)Crossover_(x,s,N_cross,Parents_num,bitmask,bitmask_num,x_dim); % 親の選択および交叉
    Mutation    = @(x,y)Mutation_(x,y,N,N_cross,N_mutat);   % 突然変異
    Selection   = @(A,Q1,Q2,S,SQ1,SQ2)Selection_(A,Q1,Q2,S,SQ1,SQ2,N,EliteCount,N_cross,N_mutat);
    J_log       = NaN(1,G+1);                               % コストログ

    % Run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 評価
    S          = Calc_Score(GroupA,N);              % 評価関数
    [J0,x0_id] = min(S);                            % スコアの最小値と要素番号
    A          = Calc_Encode(GroupA);               % GroupAをuint型に GroupA(double) -> A(uint)
    x_opt      = A(x0_id,:);                        % 最適解
    J_log(1)   = J0;
    for II = 2 : G + 1
        % 操作
        Q1  = Crossover(A,S);                       % 交叉集合
        Q2  = Mutation(A,Q1);                       % 突然変異集合
        % SQ1 = Calc_Score(Calc_Decode(Q1), N_cross); % Q1集合のスコア(1,N_cross)
        % SQ2 = Calc_Score(Calc_Decode(Q2), N_mutat); % Q2集合のスコア(1,N_mutat)
        
        % [A, fmin, fmin_id] = Selection(A,Q1,Q2,S,SQ1,SQ2,N);  % 選択
        [A, fmin, fmin_x] = Selection(A,Q1,Q2,S,[],[]);  % 選択
        S   = Calc_Score(Calc_Decode(A),N);  % 評価関数
        
        if fmin <= J0; J0 = fmin; x_opt = fmin_x; end  % 最適解の更新
        % Log
        J_log(II) = J0;
    end

    % Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = Calc_Decode(x_opt);
    LOG = J_log;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % local function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Group, fmin, fmin_x] = Selection_(A,Q1,Q2,S,SQ1,SQ2,...
                                        N,EliteCount,N_cross,N_mutat)
        % Group(A,Q1,Q2),Score(A,Q1,Q2), Group(1:N,x_dim)
        
        [fmins, Ae_id]   = mink(S,EliteCount);    % 優秀な集団を選定
        % [Q1e, Q1e_id] = mink(SQ1,EliteCount);
        % [Q2e, Q2e_id] = mink(SQ2,EliteCount);
        Ab_id = ones(1,N,"logical");
        Ab_id(Ae_id) = false;
        Ab = A(Ab_id,:);                          % エリートではない集団
        
        % ランダムに選ぶ集団
        Data_b = [Ab;Q1;Q2];                      % N - ElietCount + N_cross + N_mutat

        % 次の世代に残す集団
        Data = [A(Ae_id,:); 
            Data_b(randperm((N-ElietCount+N_cross+N_mutat),N-EliteCount),:)];     % N個選ぶ
        Group = Data(randperm(N),:);              % 偏らないようにシャッフル
        % 最適解
        fmin    = fmins(1);
        fmin_x  = A(Ae_id(1),:);
    end
    function Group = Mutation_(GroupA,GroupQ1,N,N_cross,N_mutat)
        Group = [];
    end
    function Group = Crossover_(A,Score, N_cross,Parents_num,bitmask,bitmask_num,x_dim)
        % Group(N_cross, x_dim) (N_cross = Parents_num*bitmask_num)
        % 親の選定 %
        [Dadid,Mamid] = parents_selection(Score); % Dadid(1,Parents_num)
        Dad = A(Dadid,:);                         % (Parents_num,x_dim)
        Mam = A(Mamid,:);
        % 交叉 % 
        Group = zeros(N_cross, x_dim);
        for ms_i = 1 : bitmask_num
            msk1 = bitmask(ms_i);
            msk2 = bitcmp(msk1);
            Group((Parents_num*(ms_i-1)+1):(Parents_num*ms_i),:) = ...
                    bitor(bitand(Dad(1:N_cross,:),msk1),...
                          bitand(Mam(1:N_cross,:),msk2));
        end
    end
    function [Dadid, Mamid] = Roulette_wheel(Score,Gain,Bias,N,scale,sample_num)
        % N 世代数，スケール，sample_num 親の数
        nn = Scaling(Score,bounds(Score),Gain,Bias); % (bias) < score < (gain+bias)

        % Nは無作為に抽出するためのパラメータ
        P_norm = sum(nn,"all");               % 正規化
        P      = 1 - nn./P_norm;              % 優秀なほど確率は高い
        % 親候補
        En   = fix(P*scale);
        List = NaN(1,scale+1);                % およそ(1,scale)要素
        m1   = 0;
        for i = 1 : N
            m2              = En(i) + m1;
            List(1,m1+1:m2) = i*ones(1,En(i));   % スコアが高いほど同じidがリストされる [1,1,1,2,5,5,5,5,...]
            m1              = m2;
        end
        List_num = m2;
        Boy  = randperm(List_num,sample_num); 
        Girl = randperm(List_num,sample_num);
        Dadid  = List(Boy);
        Mamid  = List(Girl);
    end

    function t = Calc_Score(Group,NN)
        t = zeros(1,NN);
        for m = 1 : NN
            t(m) = func(Group(m,:));
        end
    end
    function n = Scaling(t,xbounds,gain,bias)
        n = (t - xbounds(1))./(xbounds(2) - xbounds(1))*gain + bias; % (bias) < score < (gain+bias)
    end
    function A = Encode(GroupA,xbounds,format,Scale)
        % Group  double -> uint
        A = cast(fix(((GroupA - xbounds(1,:))./(xbounds(2,:) - xbounds(1,:)))*Scale),format);
    end
    function GroupA = Decode(A,xbounds,Scale)
        % Group  uint -> double
        GroupA = double(A)/Scale.*(xbounds(2,:) - xbounds(1,:)) + xbounds(1,:);
    end
end
%{
2024/09/30 encode,decode,完成．      calc_scoreに関しては上下限の計算を外
           に持っていくか検討中（Calc_Decode([A;Q1;Q2])の計算コストを抑え
           られるが，scoreが負となる可能性があり（親選定に障害と思ったが，
           選定時には正規化してあるので外に持っていくべき），再度正規化す
           る必要がある）
2024/10/01 Calc_Score,Roulette_wheel,Crossover_ 完成．前回のcalc_scoreはコ
　　　　　 ストと正規化を同時に行っていたので，正規化を選択方式内に移動する
           ことで計算処理を減らした．交叉に関してはペアの個数だけ子が生まれ
           ，交叉の仕方にバリエーションを加えている．したがって，交叉による
           子は1ペアから交叉の種類の数だけ生成される．Mutation_，エリート選
　　　　　 択，Selectionは全体の数が一定となるように次回設計する．
2024/10/01 Selection完成．Selection内では今後，交叉と変異のスコアも考慮した
　　　　　 選択ができるようにしてある．現状，１．もとの集団の中からエリート
           を選ぶ．２．それを除いたものと交叉と変異の集団の中からランダムに
　　　　　 選ぶ．　１と２を合わせた数が全体の個体数となるように調節した．残
　　　　　 りは，Mutation_．（交叉と変異集団のスコアを計算してもよかったが，
　　　　　 結局，次の世代で計算されるので計算しないことにした．また，スコア
　　　　　 を計算すると，局所解に陥りやすくなると考えられるため，交叉では精
　　　　　 度，変異では大域的最適性をどれだけ重要視するかが難しいと考える）

%}