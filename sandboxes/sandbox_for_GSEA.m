% sandbox_for_GSEA

%% Make Rank File
% 例のデータ
genes  = ["EphA4","Robo1","Slit2","Nrp1","Sema3a"]';
scores = [0.85, 0.67, -0.23, -0.55, 0.31]';

% NaNや欠損を除去
valid = ~isnan(scores) & genes~="";
genes  = genes(valid);
scores = scores(valid);

% （推奨）ソート：高スコア→低スコア
[~,idx] = sort(scores,'descend');
genes  = genes(idx);
scores = scores(idx);

% 書き出し（タブ区切り・ヘッダなし）
fid = fopen('ranking.rnk','w','n','UTF-8');  % UTF-8 推奨
assert(fid>0,'Cannot open file for writing.');
for i = 1:numel(genes)
    fprintf(fid,'%s\t%.10g\n', genes(i), scores(i)); % 10桁程度の有効桁
end
fclose(fid);



%%
genes = geneInfo.getGeneAcronyms();
for d = 1:dimPI
    scores = abs(relatedGeneAnalysisResults.CorrCoefMatrixSource(:,d));
    fileNameStr = strcat("geneRanking_PI",num2str(d),"_source.rnk");
    write_rnk(fileNameStr,genes,scores);
    scores = abs(relatedGeneAnalysisResults.CorrCoefMatrixTarget(:,d));
    fileNameStr = strcat("geneRanking_PI",num2str(d),"_target.rnk");
    write_rnk(fileNameStr,genes,scores);
end

%% make grp
write_grp('geneSymbols.grp',genes);


%%



%% Make Files for Normal GSEA







%% Functions
function write_rnk(filename, genes, scores, opts)
%WRITE_RNK  Write GSEA pre-ranked .rnk file (gene \t score)
%   write_rnk('out.rnk', genes, scores, opts)
%   genes: string/cellstr (N×1), scores: double (N×1)
%   opts (struct, optional fields):
%     .sort          = 'descend'|'ascend'|'' (default 'descend')
%     .dedup         = 'absmax'|'first'|'mean'|'' (default 'absmax')
%     .precision     = integer valid digits for fprintf (default 10)
%     .encoding      = 'UTF-8' (default)
%
% 例:
%   write_rnk('x.rnk', genes, scores, struct('dedup','absmax'));

if nargin < 4, opts = struct; end
if ~isfield(opts,'sort'),      opts.sort = 'descend'; end
if ~isfield(opts,'dedup'),     opts.dedup = 'absmax'; end
if ~isfield(opts,'precision'), opts.precision = 10; end
if ~isfield(opts,'encoding'),  opts.encoding = 'UTF-8'; end

genes  = string(genes(:));
scores = double(scores(:));

% remove NaN/empty
valid = ~isnan(scores) & genes~="";
genes  = genes(valid);
scores = scores(valid);

% de-duplicate
switch lower(opts.dedup)
    case 'absmax'
        [ug,~,gidx] = unique(genes, 'stable');
        keep = false(numel(genes),1);
        for k = 1:numel(ug)
            rows = find(gidx==k);
            [~,j] = max(abs(scores(rows)));
            keep(rows(j)) = true;
        end
        genes  = genes(keep); scores = scores(keep);
    case 'mean'
        [ug,~,gidx] = unique(genes, 'stable');
        s = accumarray(gidx, scores, [], @mean);
        genes = ug; scores = s;
    case 'first'
        [genes, ia] = unique(genes, 'stable');
        scores = scores(ia);
    otherwise
        % do nothing
end

% sort
switch lower(opts.sort)
    case 'descend'
        [~,idx] = sort(scores,'descend');
    case 'ascend'
        [~,idx] = sort(scores,'ascend');
    otherwise
        idx = 1:numel(scores);
end
genes  = genes(idx);
scores = scores(idx);

% write
fid = fopen(filename,'w','n',opts.encoding);
assert(fid>0, 'Cannot open %s for writing.', filename);
fmt = sprintf('%%s\\t%%.%dg\\n', opts.precision);
for i = 1:numel(genes)
    fprintf(fid, fmt, genes(i), scores(i));
end
fclose(fid);
end


function write_grp(filename, genes, varargin)
%WRITE_GRP  Write a GSEA .grp file (one gene symbol per line)
%   write_grp('set.grp', genes, 'Unique', true, 'Sort', true)
%   genes: string/cellstr vector

p = inputParser;
addParameter(p,'Unique',true,@islogical);
addParameter(p,'Sort',true,@islogical);
addParameter(p,'Encoding','UTF-8',@(x)ischar(x)||isstring(x));
parse(p,varargin{:});
opt = p.Results;

g = string(genes(:));
g = g(g~="");                 % 空を除去
if opt.Unique, g = unique(g,'stable'); end
if opt.Sort,   g = sort(g);  end

fid = fopen(filename,'w','n',opt.Encoding);
assert(fid>0,'Cannot open %s for writing.', filename);
fprintf(fid,"%s\n", g);
fclose(fid);
end