function [options,passed_on]=ScSetOptions(defaults,userargs,pass_on)
%% parses optional arguments of a function and assigns them as fields of
%% structure options
% unknown arguments are passed on into cell array passed_on if pass_on is
% present and non-empty, otherwise an error message is generated
passed_on={};
% wrap cell arguments to avoid generating multiple structs
if isstruct(defaults)
    options=defaults;
elseif iscell(defaults)
    for i=1:length(defaults)
        if iscell(defaults{i})
            defaults{i}=defaults(i);
        end
    end
    options=struct(defaults{:});
else
    error('defaults not recognized\n');
end
if nargin<3 || isempty(pass_on)
    pass_on=false;
end
for i=1:2:length(userargs)
    if isfield(options,userargs{i})
        options.(userargs{i})=userargs{i+1};
    else
        if ~pass_on
            error('option ''%s'' not recognized\n',userargs{i});
        else
            passed_on=[passed_on,{userargs{i},userargs{i+1}}]; %#ok<AGROW>
        end
    end
end
end