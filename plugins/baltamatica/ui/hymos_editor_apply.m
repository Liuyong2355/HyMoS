function path = hymos_editor_apply(text, path, save_file)
% Validate, optionally save a complete file, then commit session parameters.
check = hymos_validate_text('couette', text);
if ~check.valid
    error(['参数未应用：' char(check.message)]);
end
text = char(text);
path = char(hymos_ui_file(path));
if save_file
    path = char(hymos_ui_file(path, text));
end
had_origin = false;
try
    old_origin = evalin('base', 'hymos_config_origin');
    had_origin = true;
catch
end
origin = struct();
origin.text = text;
origin.path = path;
try
    % The text is the final commit and the sole input used by hymos_run.
    assignin('base', 'hymos_config_origin', origin);
    assignin('base', 'hymos_config_text', text);
catch err
    try
        if had_origin
            assignin('base', 'hymos_config_origin', old_origin);
        else
            evalin('base', 'clear hymos_config_origin');
        end
    catch
    end
    disp(err);
    if save_file
        error('文件已保存，但会话参数写入失败；窗口保留，请重试应用。');
    end
    error('会话参数写入失败；窗口保留，请重试应用。');
end
end
