function varargout = loadParams(params_path)
%LOADPARAMS Loads parameters from a .params file.
%
% If the file name of the .params file is not provided, the default
% image.params name is added at the end of the path.
%   Input:
%    - params_file, path to the file containing the image parameters.
%
%   Return:
%    - varargout, parameters found in the input file.
[~, ~, extension] = fileparts(params_path);

% Add the file name if not provided
if extension ~= ".params"
    params_path = params_path + "/image.params";
end

file_id = fopen(params_path);
file_content = textscan(file_id, "%s %s", 'Delimiter', '=', ...
    'CommentStyle', '%');
fclose(file_id);

varargout = cellfun(@eval, file_content{1, 2}, 'UniformOutput', false);
end