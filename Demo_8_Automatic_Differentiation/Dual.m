%{
    MATLAB上で完結する自動微分です．
    gemini 生成です．
%}
classdef Dual
    properties
        Value % 関数の値 f(x)
        Deriv % 導関数の値 f'(x)
    end
    methods
        % コンストラクタ
        function obj = Dual(value, deriv)
            if nargin == 0
                obj.Value = 0;
                obj.Deriv = 0;
            elseif nargin == 1
                obj.Value = value;
                obj.Deriv = 0; % 定数の微分は0
            else
                obj.Value = value;
                obj.Deriv = deriv;
            end
        end

        % --- 演算子のオーバーロード ---

        % 加算 (+)
        function result = plus(a, b)
            if ~isa(a, 'Dual'), a = Dual(a); end
            if ~isa(b, 'Dual'), b = Dual(b); end
            % (f+g)' = f' + g'
            result = Dual(a.Value + b.Value, a.Deriv + b.Deriv);
        end

        % 減算 (-)
        function result = minus(a, b)
            if ~isa(a, 'Dual'), a = Dual(a); end
            if ~isa(b, 'Dual'), b = Dual(b); end
            % (f-g)' = f' - g'
            result = Dual(a.Value - b.Value, a.Deriv - b.Deriv);
        end

        % 乗算 (*)
        function result = mtimes(a, b)
            if ~isa(a, 'Dual'), a = Dual(a); end
            if ~isa(b, 'Dual'), b = Dual(b); end
            % 積の微分法則: (f*g)' = f'*g + f*g'
            result = Dual(a.Value * b.Value, a.Deriv * b.Value + a.Value * b.Deriv);
        end
        function result = times(a, b)
            if ~isa(a, 'Dual'), a = Dual(a); end
            if ~isa(b, 'Dual'), b = Dual(b); end
            % 微分ルールは mtimes と同じ
            % 積の微分法則: (f*g)' = f'*g + f*g'
            result = Dual(a.Value .* b.Value, a.Deriv .* b.Value + a.Value .* b.Deriv);
            % 注意: ValueとDerivの計算を要素ごとにするため演算子を .* に変更
        end

        % 除算 (/)
        function result = mrdivide(a, b)
            if ~isa(a, 'Dual'), a = Dual(a); end
            if ~isa(b, 'Dual'), b = Dual(b); end
            % 商の微分法則: (f/g)' = (f'*g - f*g') / g^2
            result = Dual(a.Value / b.Value, (a.Deriv * b.Value - a.Value * b.Deriv) / b.Value^2);
        end
        function result = rdivide(a, b)
            if ~isa(a, 'Dual'), a = Dual(a); end
            if ~isa(b, 'Dual'), b = Dual(b); end
            % 微分ルールは mrdivide と同じ
            % 商の微分法則: (f/g)' = (f'*g - f*g') / g^2
            result = Dual(a.Value ./ b.Value, (a.Deriv .* b.Value - a.Value .* b.Deriv) ./ b.Value.^2);
            % 注意: ValueとDerivの計算を要素ごとにするため演算子を ./ と .^ に変更
        end

        % べき乗 (^)
        function result = mpower(a, b)
            if ~isa(a, 'Dual'), a = Dual(a); end
            % f(x)^c の微分は c*f(x)^(c-1)*f'(x)
            if ~isa(b, 'Dual') % bが定数の場合
                result = Dual(a.Value^b, b * a.Value^(b-1) * a.Deriv);
            else % bもDual数の場合は対数微分法などが必要で複雑になるためここでは省略
                error('Dual^Dual is not implemented.');
            end
        end

        % --- 数学関数のオーバーロード ---

        % sin
        function result = sin(a)
            % (sin(f))' = cos(f) * f'
            result = Dual(sin(a.Value), cos(a.Value) * a.Deriv);
        end

        % cos
        function result = cos(a)
            % (cos(f))' = -sin(f) * f'
            result = Dual(cos(a.Value), -sin(a.Value) * a.Deriv);
        end

        % exp
        function result = exp(a)
            % (exp(f))' = exp(f) * f'
            result = Dual(exp(a.Value), exp(a.Value) * a.Deriv);
        end

        % log
        function result = log(a)
            % (log(f))' = (1/f) * f'
            result = Dual(log(a.Value), (1/a.Value) * a.Deriv);
        end
    end
end