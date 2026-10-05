function state = hymos_editor_state(template_path)
% Read-only session restoration: unconfirmed edits never enter this state.
state = struct();
% Portable relative fallback. hymos_ui_file resolves it against the host's
% current working directory before any save; submit/run uses frozen text.
state.default_path = 'hymos_input.txt';
state.text = fileread(char(template_path));
check = hymos_validate_text('couette', state.text);
if ~check.valid
    error(['默认模板无效：' char(check.message)]);
end
state.path = state.default_path;
state.message = '已加载默认参数；应用后才生效。';
try
    saved = evalin('base', 'hymos_config_text');
catch
    return;
end
try
    check = hymos_validate_text('couette', saved);
    if ~check.valid
        state.message = '上次参数无效，已显示默认参数；原会话变量未覆盖。';
        disp(['上次参数无效：' char(check.message)]);
        return;
    end
    state.text = char(saved);
catch err
    state.message = '上次参数类型无效，已显示默认参数；原会话变量未覆盖。';
    disp(err);
    return;
end
state.message = '已恢复上次应用的参数；修改后请应用并关闭。';
try
    origin = evalin('base', 'hymos_config_origin');
    if strcmp(char(origin.text), state.text)
        state.path = char(hymos_ui_file(origin.path));
    else
        state.message = '已恢复命令窗口参数；保存位置使用当前工作目录。';
    end
catch
    % Missing, stale, protected, or unusable file metadata does not discard
    % valid applied text; use the explicit user-approved default destination.
    state.message = '已恢复参数；保存位置使用当前工作目录。';
end
end
