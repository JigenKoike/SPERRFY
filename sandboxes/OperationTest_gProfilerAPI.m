%OC_gprofilerAPI

% sandbox_gProfilerAPI

import matlab.net.*; import matlab.net.http.*;

url = URI('https://biit.cs.ut.ee/gprofiler/api/gost/profile'); % 末尾スラ無し
body_struct = struct('organism','mmusculus', ...
                     'query',{{'Gapdh'}}, ...      % 1遺伝子のみ
                     'sources',{{'GO:BP'}}, ...
                     'user_threshold',0.05);       % 余計なキーは入れない
json_body = jsonencode(body_struct);

headers = [HeaderField('Content-Type','application/json'), ...
           HeaderField('Accept','application/json'), ...
           HeaderField('User-Agent','MATLAB-gprofiler-client'), ...
           HeaderField('Expect','')];
request = RequestMessage('POST', headers, MessageBody(json_body));
resp = send(request, url);

fprintf('HTTP %s\n', string(resp.StatusCode));
raw = resp.Body.Data; if isa(raw,'uint8'), raw = char(raw.'); end
if ischar(raw), disp(raw(1:min(end,200))); end

%%
opts = weboptions('MediaType','application/json', ...
                  'Timeout', 60, ...
                  'HeaderFields', {'Accept','application/json'; ...
                                   'User-Agent','MATLAB-gprofiler-client'; ...
                                   'Expect',''});
try
    data = webwrite('https://biit.cs.ut.ee/gprofiler/api/gost/profile', body_struct, opts);
    assert(isstruct(data) && isfield(data,'result'), 'No result field.');
    fprintf('OK via webwrite; n=%d\n', numel(data.result));
catch ME
    fprintf('webwrite error: %s\n', ME.message);
end