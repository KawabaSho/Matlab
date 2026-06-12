%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Demo_8_Automatic_Differentiation
%
%   自動微分する関数を自動生成し，サンプルコードを実行します．
%
%
%   Copyright © 2025 Kawarabayashi
%   Released under the MIT license
%   https://opensource.org/licenses/mit-license.php
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 詳細な説明
%{
    　ベクトル変数x,定数a,b,cを引数とするベクトル関数f(x,a,b,c)をxで偏微分したと
    きの値を行列で出力する関数を生成します．なお，変数のサイズチェックのため，初期
    入力x0を用意してください．GenAD(@f,x0,a,b,c)を実行すると，f_jacobian.mex64
    が生成されます．ある点zにおける勾配を求める場合は，生成されたf_jacobian(z)を
    実行することで求められます．

    必須ライブラリ:
    CppAD（インストール先は'C:\Program Files (x86)\cppad'）．
    (https://cppad.readthedocs.io/stable-2025/user_guide.html）
    Matlabのデフォルトビルド環境は，Microsoft Visual C++ 2022であったので，
    Step3およびStep4はcmake --build . --target check，
    cmake --build . --target installでコマンドを通しました．cmake，grep，
    pkg-configコマンドはインストールしておいてください．
    
    必須アドオン:
    Matlab上ではMATLAB coderをアドオンしてください．
    
    生成物:
    ・codegenによるライブラリ
    ・genfuncフォルダ内に任意関数をテンプレート関数にしたソースコード
    ・自動微分する関数
    の3つです．
    　GenADで生成される自動微分関数は，引数の関数名＋_jacobian.mex64です．この
    関数の生成には，genfuncフォルダ内のコードが用いられています．生成のカスタマ
    イズはGenAD内で可能です．

    GenADアルゴリズム:
    matlab -> codegenでcppファイル生成 -> 必要箇所抽出＆ヘッダファイル生成 ->...
    mexFunction生成&内部で自動微分コード生成 -> cppでmatlab関数定義 ->...
    mexでmatlab上で動く(引数の関数名＋_jacobian)関数生成
    ※スカラー関数の自動微分のみ対応しているので，ベクトルの微分の場合は，codegen
    　ファイルで生成される形式が変化するため，mex関数に渡すcppファイルを適宜出力
    　調節する必要があります．（forを使う場合は反復回数が少ない場合は対応していま
    　すが，多くなると対応できない場合があります．）

    諸注意:
    pwdを使ってカレントディレクトリを求めているので，pwd内でスペースのあるフォルダ
    などを参照しているとエラーとなります．例 C:\source code\ は C:\source_code\
    としてください．また，勾配の流派は以下の通りです．
    f_jacobian(z)= \begin{pmatrix} \frac{\partial f_1}{\partial x_1} & \dots & \frac{\partial f_1}{\partial x_m} \\ \vdots & \ddots & \vdots \\ \frac{\partial f_n}{\partial x_1} & \dots & \frac{\partial f_n}{\partial x_m} \\ \end{pmatrix}

    引数構成：
    GenAD(@f,x0)は第一引数は関数ポインタ，第二引数は変数の初期値，第三引数以降は
    関数に与える定数値引数
    

    Example)
    GenAD(@f,x0)
    function k = f(x)
        k = 0.1*(x-1).*(x-3)+sin(x);
    end

    Example)
    GenAD(@f,x0,a)
    function k = f(x)
        k = 0.1*(x-1).*(x-3)+a*sin(x);
    end

    2025/10/22 Kawarabayahsi cssh25001@g.nihon-u.ac.jp
               (ベクトル関数のベクトル微分に対応しました．)
%}
clear,close all
x0 = 2;
GenAD(@f,x0);       % f.mの自動微分関数生成
% GenAD(@f,x0,a,b); % f(x,a,b)の引数に定数引数がある場合の書き方
%%
% 2回目の実行ではCppADの計算速度がさらに速くなります．
x0 = 2;
tic
jac = f_jacobian(x0);          % x0における自動微分実行
toc
tic
jac = 0.1*(x0-1+x0-3)+cos(x0); % f.mの導関数（理論値）
toc
tic
fk = f(Dual(x0,1));
jac = fk.Deriv;                % GeminiによるMATLABだけで完結する自動微分
toc