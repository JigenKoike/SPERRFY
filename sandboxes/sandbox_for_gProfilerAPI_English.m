%% =========================
%% 1) Input → g:Profiler request
%% =========================
% Example: gene list as cellstr (expected ~763 genes)
% gene_list = {'Dcc','Robo1','Gapdh','Vegfa'};
% When using real data (table → cell):
gene_list = table2cell(geneAcronyms_all)';

assert(iscellstr(gene_list) || isstring(gene_list), 'gene_list must be provided as cellstr or string.');
gene_list = cellstr(string(gene_list));  % normalization

body_struct = struct();
body_struct.organism = 'mmusculus';
body_struct.query = gene_list;                  % keep as is (not enclosed in curly braces)
body_struct.sources = {'GO:BP','GO:MF','GO:CC'};
body_struct.user_threshold = 0.05;              % adjust if needed
% body_struct.all_results = true;               % enable if necessary

opts = weboptions('ContentType','json','Timeout',60);
url  = 'https://biit.cs.ut.ee/gprofiler/api/gost/profile';

data = webwrite(url, body_struct, opts);
res  = data.result;
fprintf('Retrieved: %d (GO terms)\n', numel(res));
assert(~isempty(res), 'Result is empty. Please review parameters.');


%% =========================
%% 2) Matching 763 ↔ 739 genes (alignment and missing detection)
%% =========================
% 2-1) Expand the result table vertically (convert 1-row table with {739x1 cell} columns → 739-row table)
T = struct2table(res,'AsArray',true);   % often height=1
if height(T) == 1
    vars = T.Properties.VariableNames;
    % First, estimate N (= number of recognized genes) from intersections
    assert(iscell(T.intersections) && isscalar(T.intersections{1}), 'Unexpected shape of intersections column.');
    N = numel(T.intersections{1});
    % Empty container
    T2 = table('Size',[N numel(vars)], ...
               'VariableTypes', repmat("cell",1,numel(vars)), ...
               'VariableNames', vars);
    % Expand each column
    for i = 1:numel(vars)
        c = T.(vars{i});
        if iscell(c) && isscalar(c)
            inner = c{1};
            if isrow(inner), inner = inner.'; end
            if size(inner,1) == N
                T2.(vars{i}) = inner;
            else
                % Broadcast scalar etc. to all rows
                T2.(vars{i}) = repmat(c(1), N, 1);
            end
        else
            T2.(vars{i}) = repmat(c(1), N, 1);
        end
    end
else
    T2 = T;  % already expanded vertically
end
fprintf('Expanded table size: %d x %d\n', height(T2), width(T2));
assert(ismember("intersections", string(T2.Properties.VariableNames)), 'intersections column not found.');

%% 2-2) Retrieve “recognized 739 genes” from meta (Ensembl → Symbol)
gm = data.meta.genes_metadata;
ensgs = string(gm.query.query_1.ensgs);     % 739x1 Ensembl IDs (same order as intersections)
mapS  = gm.query.query_1.mapping;           % struct: field name = Symbol (contains Ensembl)
failed_list = string(gm.failed);            % officially unrecognized inputs (e.g., 22 genes)

% Create reverse lookup map Ensembl → Symbol
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
fprintf('Reverse Ensembl→Symbol mapping: %d entries\n', ensg2sym.Count);

% Symbol array in 739 order (recognized_syms)
recognized_syms = strings(numel(ensgs),1);
for i = 1:numel(ensgs)
    k = char(ensgs(i));
    if isKey(ensg2sym,k), recognized_syms(i) = string(ensg2sym(k)); else, recognized_syms(i) = ""; end
end

%% 2-3) Match against the 763 input (positions and missing)
gene_input = string(gene_list(:));
GI = upper(gene_input);
GR = upper(recognized_syms);

[tf_rec_in_inp, pos_in_input] = ismember(GR, GI);     % 739 → 763 placement
[tf_inp_in_rec, ~] = ismember(GI, GR);                % which of 763 were recognized
missing_genes = gene_input(~tf_inp_in_rec);           % unrecognized genes (estimated)
fprintf('Input %d / Recognized %d / Unrecognized (est.) %d / API failed %d\n', ...
        numel(gene_input), numel(recognized_syms), numel(missing_genes), numel(failed_list));

%% =========================
%% 3) Create gene × GO logical matrix
%% =========================
% 3-1) First make 739×GO (intersections: for each GO row, a 739×1 cell)
nTerms = height(T2);
nGenes739 = numel(T2.intersections{1});   % should be 739
hit_any = false(nGenes739, nTerms);       % non-empty = hit
hit_exp = false(nGenes739, nTerms);       % hit only if includes experimental evidence (IMP/IDA/IPI/IGI/IEP/EXP)
EXP = ["IMP","IDA","IPI","IGI","IEP","EXP"];

for t = 1:nTerms
    evid_col = T2.intersections{t};      % 739×1 cell ([], 'IEA', {'IEA','IMP'}, etc.)
    if isrow(evid_col), evid_col = evid_col.'; end
    m = min(nGenes739, numel(evid_col));
    for g = 1:m
        x = evid_col{g};
        % non-empty means hit
        hit_any(g,t) = ~isempty(x);
        % whether it contains experimental evidence codes (handle nested)
        if ~isempty(x)
            % normalize x into array of string codes (inline processing)
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
fprintf('Created: hit_any=%dx%d, hit_exp=%dx%d\n', size(hit_any), size(hit_exp));

%% 3-2) Convert 739×GO → 763×GO (restore original gene_list order)
[nRec, nTerms] = size(hit_any);
assert(nRec == numel(recognized_syms), 'Row count mismatch between 739 matrix and recognized_syms.');

hit_any_full = false(numel(gene_input), nTerms);
hit_exp_full = false(numel(gene_input), nTerms);
valid = tf_rec_in_inp & (pos_in_input > 0);
hit_any_full( pos_in_input(valid), : ) = hit_any( valid, : );
hit_exp_full( pos_in_input(valid), : ) = hit_exp( valid, : );

fprintf('All-zero rows (any evidence): %d, all-zero rows (experimental only): %d\n', ...
        sum(~any(hit_any_full,2)), sum(~any(hit_exp_full,2)));

%% 3-3) Sanitize column names (avoid overly long ones)
% 1) Use short IDs (term_id/native) as column names
if ismember("term_id", string(T2.Properties.VariableNames))
    termIDs = string(T2.term_id);
elseif ismember("native", string(T2.Properties.VariableNames))
    termIDs = string(T2.native);
else
    % If no ID is available, generate short placeholder names (e.g., GO_0001)
    termIDs = "GO_" + string(1:nTerms);
end

% 2) Full names (term_name/name) go to VariableDescriptions
if ismember("term_name", string(T2.Properties.VariableNames))
    termNamesFull = string(T2.term_name);
elseif ismember("name", string(T2.Properties.VariableNames))
    termNamesFull = string(T2.name);
else
    termNamesFull = strings(nTerms,1);
end

% 3) Make valid MATLAB variable names (handle length and duplicates)
safeVarNames = cellstr( matlab.lang.makeValidName(termIDs) );
safeVarNames = matlab.lang.makeUniqueStrings(safeVarNames, {}, namelengthmax);

% 4) Create tables (column name = sanitized ID, row name = original gene_list)
HitAnyFull = array2table(hit_any_full, 'RowNames', cellstr(gene_input), 'VariableNames', safeVarNames);
HitExpFull = array2table(hit_exp_full, 'RowNames', cellstr(gene_input), 'VariableNames', safeVarNames);

% 5) Attach full term names as descriptions
HitAnyFull.Properties.VariableDescriptions = cellstr(termNamesFull);
HitExpFull.Properties.VariableDescriptions = cellstr(termNamesFull);

% 6) Also prepare correspondence table (column name ↔ display name) for readability
GO_legend = table(safeVarNames(:), termIDs(:), termNamesFull(:), ...
    'VariableNames', {'VarName','TermID','TermName'});

% Display preview
disp(HitAnyFull(1:min(5,height(HitAnyFull)), 1:min(5,width(HitAnyFull))));
disp(GO_legend(1:min(5,height(GO_legend)), :));

%% =========================
%% 4) Extract specific GO(s) → logical array (optional: limit to BP, etc.)
%% =========================
% Input example:
%targets = ["axon guidance","metabolic","angiogenesis","GO:0006955"];  % mix of names and IDs allowed
targets = ["synap"];
matchMode = "contains";     % "contains" or "exact" (applied to names)
limitSource = false;        % true to restrict to specific GO categories (e.g., GO:BP)
wantedSources = "GO:BP";    % e.g., only biological process

% --- If legend lacks Source column, add from T2 ---
if exist('T2','var') && ~ismember("Source", string(GO_legend.Properties.VariableNames))
    GO_legend.Source = string(T2.source);
end

% --- Prepare search keys ---
names = string(GO_legend.TermName);
ids   = string(GO_legend.TermID);
vars  = string(GO_legend.VarName);

% --- Determine which columns to extract (ID prioritized; fall back to name) ---
take = false(height(GO_legend),1);
for q = targets(:).'
    s = string(q);
    if startsWith(upper(s),"GO:")
        % ID-based (exact match only)
        take = take | (upper(ids) == upper(s));
    else
        % Name-based
        switch lower(matchMode)
            case "exact"
                take = take | (upper(names) == upper(s));
            otherwise % contains
                take = take | contains(upper(names), upper(s));
        end
    end
end

% --- Filter by source if specified ---
if limitSource && ismember("Source", string(GO_legend.Properties.VariableNames))
    take = take & ismember(string(GO_legend.Source), string(wantedSources));
end

colIdx = find(take);
if isempty(colIdx)
    error('No columns matched the specified GO terms. Please review targets.');
end

% --- Extract selected columns from HitAnyFull / HitExpFull ---
selVars = cellstr(vars(colIdx));               % valid column names (safe)
HitsAny_sel = HitAnyFull(:, selVars);
HitsExp_sel = HitExpFull(:, selVars);

% --- Also obtain logical arrays (rows = genes, cols = selected GO terms) ---
L_any = table2array(HitsAny_sel);              % logical (Ngenes x Ncols)
L_exp = table2array(HitsExp_sel);              % logical (Ngenes x Ncols)

% --- Column descriptions (readable term names or IDs) ---
Legend_sel = GO_legend(colIdx, :);
% Readable output (preview of first few)
fprintf('Extracted columns: %d (rows=%d genes)\n', numel(colIdx), size(L_any,1));
disp(Legend_sel(1:min(8,height(Legend_sel)), :));

% Example: extract 0/1 vector of the first GO column
v_any_first = L_any(:,1);  % logical vector in same order as gene_list