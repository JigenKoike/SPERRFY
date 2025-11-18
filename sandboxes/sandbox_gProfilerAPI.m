%% =========================
%% 1) 入力 → g:Profiler 呼び出し
%% =========================
% 例: cellstr の遺伝子リスト（763件想定）
% gene_list = {'Dcc','Robo1','Gapdh','Vegfa'};
% 実データを使うとき（テーブル→cell）:
geneAcronyms_all = geneInformationTable(:,'acronym');
gene_list = table2cell(geneAcronyms_all)';

assert(iscellstr(gene_list) || isstring(gene_list), 'gene_list は cellstr/string で与えてください。');
gene_list = cellstr(string(gene_list));  % 正規化

body_struct = struct();
body_struct.organism = 'mmusculus';
body_struct.query = gene_list;                  % ← そのまま（波括弧で包まない）
body_struct.sources = {'GO:BP','GO:MF','GO:CC'};
body_struct.user_threshold = 0.05;              % 必要に応じて調整
% body_struct.all_results = true;               % 必要なら有効化

opts = weboptions('ContentType','json','Timeout',60);
url  = 'https://biit.cs.ut.ee/gprofiler/api/gost/profile';

data = webwrite(url, body_struct, opts);
res  = data.result;
fprintf('取得件数: %d（GO terms）\n', numel(res));
assert(~isempty(res), '結果が空です。パラメータを見直してください。');


%% =========================
%% 2) 763 ↔ 739 の対応づけ（順序の確定と欠落の把握）
%% =========================
% 2-1) result を行方向に展開（各列が {739x1 cell} を抱えている 1行テーブル → 739行テーブルへ）
T = struct2table(res,'AsArray',true);   % 高さ1のことが多い
if height(T) == 1
    vars = T.Properties.VariableNames;
    % まず N（= 認識遺伝子数）を intersections から推定
    assert(iscell(T.intersections) && isscalar(T.intersections{1}), 'intersections 列の形が想定外です。');
    N = numel(T.intersections{1});
    % 空の器
    T2 = table('Size',[N numel(vars)], ...
               'VariableTypes', repmat("cell",1,numel(vars)), ...
               'VariableNames', vars);
    % 各列を展開
    for i = 1:numel(vars)
        c = T.(vars{i});
        if iscell(c) && isscalar(c)
            inner = c{1};
            if isrow(inner), inner = inner.'; end
            if size(inner,1) == N
                T2.(vars{i}) = inner;
            else
                % スカラー等は全行にブロードキャスト
                T2.(vars{i}) = repmat(c(1), N, 1);
            end
        else
            T2.(vars{i}) = repmat(c(1), N, 1);
        end
    end
else
    T2 = T;  % すでに行方向ならそのまま
end
fprintf('展開後のサイズ: %d x %d\n', height(T2), width(T2));
assert(ismember("intersections", string(T2.Properties.VariableNames)), 'intersections 列が見つかりません。');

%% 2-2) g:Profiler が「認識した 739 遺伝子の並び」を meta から取得（Ensembl → Symbol）
gm = data.meta.genes_metadata;
ensgs = string(gm.query.query_1.ensgs);     % 739x1 Ensembl ID（intersections の並びと一致）
mapS  = gm.query.query_1.mapping;           % struct: フィールド名=Symbol（中に Ensembl が入っている）
failed_list = string(gm.failed);            % 公式に未認識とされた入力（例: 22件）

% 逆引きマップ Ensembl → Symbol を作成
ensg2sym = containers.Map('KeyType','char','ValueType','char');
fns = fieldnames(mapS);
for i = 1:numel(fns)
    sym = fns{i}; entry = mapS.(sym);
    if iscell(entry), entry = entry{1}; end
    ensg = '';
    if isstruct(entry)
        subf = fieldnames(entry);
        for k = 1:numel(subf)
            val = entry.(subf{k});
            if ischar(val) || (isstring(val) && isscalar(val))
                s = char(string(val));
                if ~isempty(regexp(s,'^ENS\w+G\d+$','once')), ensg = s; break; end
            end
        end
    elseif ischar(entry) || (isstring(entry) && isscalar(entry))
        s = char(string(entry));
        if ~isempty(regexp(s,'^ENS\w+G\d+$','once')), ensg = s; end
    end
    if ~isempty(ensg), ensg2sym(ensg) = sym; end
end
fprintf('Ensembl→Symbol 逆引き: %d 件\n', ensg2sym.Count);

% 739 並びのシンボル名配列（recognized_syms）
recognized_syms = strings(numel(ensgs),1);
for i = 1:numel(ensgs)
    k = char(ensgs(i));
    if isKey(ensg2sym,k), recognized_syms(i) = string(ensg2sym(k)); else, recognized_syms(i) = ""; end
end

%% 2-3) 763 入力と突き合わせ（対応位置・欠落）
gene_input = string(gene_list(:));
GI = upper(gene_input);
GR = upper(recognized_syms);

[tf_rec_in_inp, pos_in_input] = ismember(GR, GI);     % 739 → 763 の置き場所
[tf_inp_in_rec, ~] = ismember(GI, GR);                % 763 のうち認識されたもの
missing_genes = gene_input(~tf_inp_in_rec);           % 認識されなかった遺伝子（推定）
fprintf('入力 %d / 認識 %d / 未認識(推定) %d / API failed %d\n', ...
        numel(gene_input), numel(recognized_syms), numel(missing_genes), numel(failed_list));

%% =========================
%% 3) 遺伝子×GO の logical 行列を作成
%% =========================
% 3-1) まず 739×GO を作る（intersections: 各 GO 行に 739×1 cell）
nTerms = height(T2);
nGenes739 = numel(T2.intersections{1});   % = 739 のはず
hit_any = false(nGenes739, nTerms);       % 非空=ヒット
hit_exp = false(nGenes739, nTerms);       % 実験系（IMP/IDA/IPI/IGI/IEP/EXP）を含む時のみヒット
EXP = ["IMP","IDA","IPI","IGI","IEP","EXP"];

for t = 1:nTerms
    evid_col = T2.intersections{t};      % 739×1 cell（[], 'IEA', {'IEA','IMP'} など）
    if isrow(evid_col), evid_col = evid_col.'; end
    m = min(nGenes739, numel(evid_col));
    for g = 1:m
        x = evid_col{g};
        % 非空ならヒット
        hit_any(g,t) = ~isempty(x);
        % 実験系コードを含むか（入れ子に対応）
        if ~isempty(x)
            % x を「文字列コードの配列」に正規化（関数なしでその場処理）
            codes = strings(0,1);
            stk = {x};
            while ~isempty(stk)
                y = stk{1}; stk(1) = [];
                if iscell(y)
                    stk = [y(:).', stk]; %#ok<AGROW>
                elseif ischar(y) || (isstring(y) && isscalar(y))
                    codes(end+1,1) = string(y); %#ok<AGROW>
                end
            end
            hit_exp(g,t) = any(ismember(upper(codes), EXP));
        end
    end
end
fprintf('作成: hit_any=%dx%d, hit_exp=%dx%d\n', size(hit_any), size(hit_exp));

%% 3-2) 739×GO → 763×GO へ復元（元の gene_list 順）
[nRec, nTerms] = size(hit_any);
assert(nRec == numel(recognized_syms), '739行列と recognized_syms の行数が一致しません。');

hit_any_full = false(numel(gene_input), nTerms);
hit_exp_full = false(numel(gene_input), nTerms);
valid = tf_rec_in_inp & (pos_in_input > 0);
hit_any_full( pos_in_input(valid), : ) = hit_any( valid, : );
hit_exp_full( pos_in_input(valid), : ) = hit_exp( valid, : );

fprintf('全0行（非空=ヒット基準）: %d, 全0行（実験系のみ）: %d\n', ...
        sum(~any(hit_any_full,2)), sum(~any(hit_exp_full,2)));

%% 3-3) 列名を安全化（長い名前は不可）
% 1) まず短いID候補（term_id/native）を列名に使う
if ismember("term_id", string(T2.Properties.VariableNames))
    termIDs = string(T2.term_id);
elseif ismember("native", string(T2.Properties.VariableNames))
    termIDs = string(T2.native);
else
    % IDがない場合は短縮名を生成（GO_0001 のように）
    termIDs = "GO_" + string(1:nTerms);
end

% 2) 説明用のフル名（term_name/name）を VariableDescriptions へ
if ismember("term_name", string(T2.Properties.VariableNames))
    termNamesFull = string(T2.term_name);
elseif ismember("name", string(T2.Properties.VariableNames))
    termNamesFull = string(T2.name);
else
    termNamesFull = strings(nTerms,1);
end

% 3) VariableNames 用に MATLAB 有効名へ整形（長さ・重複にも対応）
safeVarNames = cellstr( matlab.lang.makeValidName(termIDs) );
safeVarNames = matlab.lang.makeUniqueStrings(safeVarNames, {}, namelengthmax);

% 4) 表を作成（列名=安全化ID、行名=元の gene_list）
HitAnyFull = array2table(hit_any_full, 'RowNames', cellstr(gene_input), 'VariableNames', safeVarNames);
HitExpFull = array2table(hit_exp_full, 'RowNames', cellstr(gene_input), 'VariableNames', safeVarNames);

% 5) 説明（フル term 名）を列に付与
HitAnyFull.Properties.VariableDescriptions = cellstr(termNamesFull);
HitExpFull.Properties.VariableDescriptions = cellstr(termNamesFull);

% 6) 便利用の対応表（列名↔表示名）も作っておくと後で見やすい
GO_legend = table(safeVarNames(:), termIDs(:), termNamesFull(:), ...
    'VariableNames', {'VarName','TermID','TermName'});

% 確認表示
disp(HitAnyFull(1:min(5,height(HitAnyFull)), 1:min(5,width(HitAnyFull))));
disp(GO_legend(1:min(5,height(GO_legend)), :));

%% =========================
%% 4) 指定GOだけ抽出 → logical配列（任意: BPに限定など）
%% =========================
% 入力（例）:
%targets = ["axon guidance","metabolic","angiogenesis","GO:0006955"];  % 名前とID混在OK
targets = ["synap"];
matchMode = "contains";     % "contains" or "exact"（名前に対して）
limitSource = false;        % true にすると GO:BP などに限定
wantedSources = "GO:BP";    % 例: 生物学的過程のみ

% --- legend に Source 列が無ければ T2 から付ける ---
if exist('T2','var') && ~ismember("Source", string(GO_legend.Properties.VariableNames))
    GO_legend.Source = string(T2.source);
end

% --- 探索キーを用意 ---
names = string(GO_legend.TermName);
ids   = string(GO_legend.TermID);
vars  = string(GO_legend.VarName);

% --- どの列を取るか決める（ID優先。ID指定が無ければ名前で検索） ---
take = false(height(GO_legend),1);
for q = targets(:).'
    s = string(q);
    if startsWith(upper(s),"GO:")
        % ID指定（厳密一致のみ）
        take = take | (upper(ids) == upper(s));
    else
        % 名前指定
        switch lower(matchMode)
            case "exact"
                take = take | (upper(names) == upper(s));
            otherwise % contains
                take = take | contains(upper(names), upper(s));
        end
    end
end

% --- ソースで絞る（任意） ---
if limitSource && ismember("Source", string(GO_legend.Properties.VariableNames))
    take = take & ismember(string(GO_legend.Source), string(wantedSources));
end

colIdx = find(take);
if isempty(colIdx)
    error('指定した GO に一致する列が見つかりませんでした。targets を見直してください。');
end

% --- HitAnyFull / HitExpFull から該当列だけ抜き出す ---
selVars = cellstr(vars(colIdx));               % 有効列名（safe）
HitsAny_sel = HitAnyFull(:, selVars);
HitsExp_sel = HitExpFull(:, selVars);

% --- logical 配列としても取得（行=gene、列=選ばれたGO） ---
L_any = table2array(HitsAny_sel);              % logical (Ngenes x Ncols)
L_exp = table2array(HitsExp_sel);              % logical (Ngenes x Ncols)

% --- 列の説明（人が読む名称やID） ---
Legend_sel = GO_legend(colIdx, :);
% 可読出力（先頭だけプレビュー）
fprintf('抽出列: %d 本（行=%d genes）\n', numel(colIdx), size(L_any,1));
disp(Legend_sel(1:min(8,height(Legend_sel)), :));

% 例）個別ベクトルの取り出し：最初のGO列の 0/1 ベクトル
v_any_first = L_any(:,1);  % gene_list と同じ順の logical