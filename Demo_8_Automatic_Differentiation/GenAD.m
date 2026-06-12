%{
    自動微分する関数を自動生成できます．詳しい使い方はMain.mをご確認ください．

    依存関係：
    1. CppAD（インストール先は'C:\Program Files (x86)\cppad'）．
      (https://cppad.readthedocs.io/stable-2025/user_guide.html）
    2. MATLAB coder
    
    ex)
        ベクトル変数x,定数a,b,cのベクトル関数f(x,a,b,c)をxで偏微分したときの値を
    行列で出力する関数を生成します．なお，変数のサイズチェックのため，初期入力x0を
    用意してください．GenAD(@f,x0,a,b,c)を実行すると，f_jacobian.mex64が生成され
    ます．ある点zにおける勾配を求める場合は，f_jacobian(z)で求められます．勾配を求
    める際に，定数引数を使いたい場合は，
        line 47: coder.Constant，
        convert_to_template_nのtemplateCode，
        generate_template_cpp
    を変更する必要があります．

    勾配の流派は以下の通りです．
    f_jacobian(z)= \begin{pmatrix} \frac{\partial f_1}{\partial x_1} & \dots & \frac{\partial f_1}{\partial x_m} \\ \vdots & \ddots & \vdots \\ \frac{\partial f_n}{\partial x_1} & \dots & \frac{\partial f_n}{\partial x_m} \\ \end{pmatrix}

  Copyright © 2025 kawarabayashi
  Released under the MIT license
  https://opensource.org/licenses/mit-license.php
%}
function GenAD(infunc,x0,varargin)
    % x0 はinfuncの引数の変数としての入力する値
    % vararginは定数値引数の引き渡し
    %{
      <infunc内のルール>
        ・^2はサポートできますが，^3はサポート外です．
        ・うまくいかない場合は，codegen/lib/infunc.cppを確認してどのようなものが
        　生成され，genfunc/infunc.hとinfunc_jacobian.cppでどのように生成されて
          いるか比較すること．同値であれば動作します．
        ・infuncの出力はスカラー量でお願いします．今後，対応させます．
          -> 対応しました．2025/10/22 kawarabayashi cssh25001@g.nihon-u.ac.jp
    %}
    cfg = coder.config('lib');
    cfg.TargetLang = 'C++';
    cfg.InstructionSetExtensions = 'None';
    functionName = func2str(infunc);
    m = length(x0);
    if isempty(varargin) % nargin
        n = length(infunc(x0));
        codegen(functionName,'-config',cfg,'-args',{x0})
    else
        n = length(infunc(x0,varargin{:}));
        navaragin = nargin - 2;
        in_x = cell(1,navaragin+1);
        in_x{1} = x0;
        for i = 1 : navaragin
            in_x{i+1} = coder.Constant(varargin{i});
        end
        codegen(functionName,'-config',cfg,'-args',in_x)
    end
    % 以上で@(x)f(x,constant values)をC++で生成

    % template関数として再定義
    gendir        = "genfunc";
    path_genfile  = pwd + "\" + gendir;
    genid         = pwd + "\codegen\lib\" + functionName + "\" + functionName + ".cpp";
    genout        = path_genfile + "\" + functionName + ".hpp";
    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    mkdir(gendir);
    warning('on', 'MATLAB:MKDIR:DirectoryExists');
    
    % genfuncフォルダに.hppファイルを生成
    generate_code(genout,convert_to_template_n(genid, functionName,m,n))
    fprintf("\t%s.hpp was generated.\n",functionName)

    % cppadの組み込み
    outgen  = path_genfile + "\" + functionName+"_jacobian.cpp"; % cppadをmatlabで呼び出すためのcpp生成ファイルのディレクトリ
    IncHeader = gendir + "/" + functionName + ".hpp";

    generate_code(outgen, generate_template_cpp, IncHeader, num2str(m), num2str(n), functionName);
    fprintf("\t%s_jacobian.cpp was generated.\n\n",functionName)
    cppad_include_path = 'C:\Program Files (x86)\cppad\include'; % インストールしたcppadのヘッダーファイル読み込み
    mex(['-I' cppad_include_path],['-I' convertStringsToChars(pwd)], outgen)
end
% .hpp for cppad
function templateCode = convert_to_template_n(generatedFile, templateFuncName,m,n)
    % generatedFile: codegenが生成したファイル名 (例: 'f.cpp')
    % templateFile:  テンプレートとして保存する新しいファイル名 (例: 'my_func.h')
    % templateFuncName: 新しいテンプレート関数名 (例: 'my_func')

    % 1. 生成されたC++ファイルを読み込む
    code = fileread(generatedFile);

    % 予約関数の保護
    code = ProtectFunction(code); % std::cosをcosに変換したり，変数x以外の文字を別の文字に一度置き換えたりしている
    
    % 2. 正規表現を使って関数定義部分を抽出・置換する
    % 'real_T f(real_T x)' のような形式を検索
    p0   = strfind(code,"// Include Files");
    p    = strfind(code,"{");
    pend = strfind(code,"}");
    if n == 1 
        % infuncがスカラー量のときは，codegenでreturnを用いてcppファイルが生成されるため
        % cppadで必要な配列での形式（例 output[0]）に変換する必要がある．
        outputArrayName = "y_func_mask01i";
        p_return    = strfind(code,"return");
        p_mainstart = p_return(1) + 6;           % 関数内部の抽出のためのポインタ
        ifOutIsScalar = "y[0] =";                % output[0]
    else
        p_arg0 = strfind(code,"//                double") + strlength("//                double")+1; % 出力の名前の先頭
        p_arg1 = strfind(code(p_arg0:end),"[")+strlength(code(1:p_arg0))-2; % 出力の名前の最後尾
        outputArrayName = code(p_arg0(1):p_arg1(1));
        p_mainstart = p(1)+3;
        ifOutIsScalar = "";
    end

    convert_code_xarray = code(p_mainstart:pend(end)-2); % 関数内部
    convert_code_xarray = strrep(convert_code_xarray, outputArrayName, 'funcmask0p0z');
    % スカラー引数の場合，codegenでは配列として生成されないので，ここで配列に書き換える必要がある
    if m == 1
        convert_code_xarray = strrep(convert_code_xarray, 'x', 'x[0]');
    end

    % 予約関数の適用
    if n == 1
        outputArrayName = "y";
    end
    convert_code_xarray = ApplyFunction(convert_code_xarray);
    convert_code_xarray = strrep(convert_code_xarray, 'funcmask0p0z', outputArrayName);

    % 3. テンプレート形式のコードを生成
    templateCode = [
        code(1:p0(1)-5)
        ""
        "#pragma once"
        "#include <cmath>"
        "#include <cppad/cppad.hpp>"
        ""
        "template <typename T> "+ "void" + " " + ...
        templateFuncName + "(const CppAD::vector<T>& " + "x" + ", " +...
        "CppAD::vector<T>& " + outputArrayName + ")"
        "{"
        ifOutIsScalar + convert_code_xarray
        "}"
    ];
    templateCode = strjoin(templateCode, "\n");
end
function code = ProtectFunction(code)
    % 関数本体内の 'muDoubleScalarSin' などをcppadに対応できる 'sin' に置換
    code = strrep(code, 'std::sin', 'sin');
    code = strrep(code, 'std::cos', 'cos');
    code = strrep(code, 'std::pow', 'pow');
    code = strrep(code, 'std::sqrt', 'sqrt');
    code = strrep(code, 'std::log',  'log');
    code = strrep(code, 'std::exp', 'std::eep');%mask 入力がスカラーのとき，配列として書き直す必要があるので，cmathのxが付く関数をマスキングする．
    % 他にもあればここに追加
end
function code = ApplyFunction(code)
    code = strrep(code, 'std::eep', 'exp');% mask
    code = strrep(code, 'double', 'T');
    code = strrep(code, 'int', 'size_t');
end
function generate_code(templateFile,templateCode,varargin)
    % 4. 書き出す
    % templateFile 出力先ファイル, templateCode 保存したい文字列, varargin　変数文字列の挿入
    fid = fopen(templateFile, 'w');
    fprintf(fid, templateCode, varargin{:});
    fclose(fid);
end

% jacobian.cpp for matlab
function cpp_code_str = generate_template_cpp()

tmp_cpp_code = [
"// func_jacobian.cpp"
"#include ""mex.h"""
"#include <cppad/cppad.hpp>"
"#include <cmath>"
"#include ""%s"""                                                          % codegenで生成したヘッダーファイル
""
""
"void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])"
"{"
"    "
"    double *x_in = mxGetPr(prhs[0]);"
"    "
"    using CppAD::AD;"
"    size_t n = %s;"                                                       % 入力次元
"    size_t m = %s;"                                                       % ヤコビアンの出力次元
"    CppAD::vector<AD<double>> X(n);"
"    for(size_t i = 0; i < n; i++) {"
"        X[i] = x_in[i];"
"     }"
"    "
"    CppAD::Independent(X);"
"    "
"    CppAD::vector<AD<double>> Y(m);"
"    %s(X,Y);"                                                             % function name 
"    "
"    CppAD::ADFun<double> F(X,Y);" 
"    "
"    CppAD::vector<double> x(n);"
"    for(size_t i = 0; i < n; i++) {"
"       x[i] = x_in[i];"
"    }"
"    "
"    CppAD::vector<double> jacobian = F.Jacobian(x);"
"    "
"    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);"
"    double *y_out = mxGetPr(plhs[0]);"
"    "
"    const double *jac_ptr = &jacobian[0];"
"    for(size_t j = 0; j < n; j++) {"
"        for(size_t i = 0; i < m; i++) {"
"            size_t dest_index = j * m + i;"
"            size_t src_index = i * n + j;"
"            y_out[dest_index] = jac_ptr[src_index];"
"        }"
"     }"
"}"
];

cpp_code_str = strjoin(tmp_cpp_code, "\n");
end
