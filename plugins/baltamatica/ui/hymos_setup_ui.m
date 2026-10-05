function fig = hymos_setup_ui(case_id, template_path)
% A plain category form backed by the complete original TXT/INI draft.
if ~strcmpi(char(case_id), 'couette')
    error('hymos_setup currently supports only couette');
end
state = hymos_editor_state(template_path);
ui.schema = hymos_form_schema();
fig = uifigure('Name', 'HyMoS 一维平行板 BGK 参数配置', 'Position', [60 45 1180 820]);
% Light sidebar, large content title, and a compact fixed action bar.
hymos_visual_style(fig, 'Color', [0.97 0.98 1]);
% This UI uses a fixed pixel layout; prevent a smaller window from making
% labels, hints, and input fields overlap on hosts without responsive grids.
hymos_visual_style(fig, 'Resize', 'off');
sidebar = uilabel(fig, 'Text', '', 'Position', [0 155 220 665]);
hymos_visual_style(sidebar, 'BackgroundColor', [0.94 0.96 0.99]);
brand = uilabel(fig, 'Text', 'HyMoS', 'Position', [28 753 175 38]);
hymos_visual_style(brand, 'FontSize', 25, 'FontWeight', 'bold', 'FontColor', [0.09 0.22 0.43]);
brand_note = uilabel(fig, 'Text', '1D PARALLEL PLATE / BGK', 'Position', [28 720 180 25]);
hymos_visual_style(brand_note, 'FontSize', 13, 'FontColor', [0.39 0.46 0.58]);
ui.title = uilabel(fig, 'Text', '一维平行板 BGK 算例', 'Position', [260 746 875 45]);
hymos_visual_style(ui.title, 'FontSize', 30, 'FontWeight', 'bold', 'FontColor', [0.08 0.16 0.30]);
ui.note = uilabel(fig, 'Text', '选择左侧分类，配置本次计算参数。', 'Position', [262 704 875 30]);
hymos_visual_style(ui.note, 'FontSize', 16, 'FontColor', [0.40 0.46 0.56]);
ui.section = uilabel(fig, 'Text', '', 'Position', [260 650 875 42]);
hymos_visual_style(ui.section, 'FontSize', 20, 'FontWeight', 'bold', ...
                   'FontColor', [0.12 0.30 0.58], 'BackgroundColor', [0.91 0.95 1]);
ui.editor = uitextarea(fig, 'Position', [260 180 875 462], ...
                      'Value', hymos_form_lines(state.text));
hymos_visual_style(ui.editor, 'FontSize', 17);
ui.path = uilabel(fig, 'Text', ['当前文件：' state.path], 'Position', [25 117 1115 25]);
hymos_visual_style(ui.path, 'FontSize', 14, 'FontColor', [0.30 0.37 0.48]);
ui.message = uilabel(fig, 'Text', state.message, 'Position', [260 153 875 25]);
hymos_visual_style(ui.message, 'FontSize', 14, 'FontColor', [0.30 0.39 0.53]);
ui.save_on_apply = uicheckbox(fig, 'Text', '应用时保存到当前文件（直接覆盖）', ...
                            'Value', false, 'Position', [25 79 570 28]);
hymos_visual_style(ui.save_on_apply, 'FontSize', 15);
ui.nav = cell(1, 5);
for i = 1:5
    ui.nav{i} = uibutton(fig, 'Text', ui.schema.groups{i}, ...
                         'Position', [18 625-(i-1)*76 185 56]);
    hymos_visual_style(ui.nav{i}, 'FontSize', 19, 'FontColor', [0.28 0.34 0.44], ...
                       'BackgroundColor', [0.94 0.96 0.99]);
end
load_button = uibutton(fig, 'Text', '加载参数', 'Position', [25 20 130 45]);
save_button = uibutton(fig, 'Text', '保存参数', 'Position', [170 20 130 45]);
save_as_button = uibutton(fig, 'Text', '另存为', 'Position', [315 20 115 45]);
default_button = uibutton(fig, 'Text', '恢复默认', 'Position', [445 20 130 45]);
validate_button = uibutton(fig, 'Text', '检查参数', 'Position', [590 20 140 45]);
apply_button = uibutton(fig, 'Text', '应用并关闭', 'Position', [925 20 210 45]);
buttons = {load_button, save_button, save_as_button, default_button, validate_button};
for i = 1:length(buttons)
    hymos_visual_style(buttons{i}, 'FontSize', 17, 'BackgroundColor', [1 1 1], ...
                       'FontColor', [0.17 0.25 0.38]);
end
hymos_visual_style(apply_button, 'FontSize', 18, 'FontWeight', 'bold', ...
                   'BackgroundColor', [0.10 0.36 0.88], 'FontColor', [1 1 1]);
ui.labels = cell(1, 8); ui.edits = cell(1, 8);
ui.choices = cell(1, 8); ui.hints = cell(1, 8); ui.extra = {};
ui.advanced_help = [];
ui.form_available = true;
try
    ui.type = uilabel(fig, 'Text', '当前工况：Couette', 'Position', [280 600 830 36]);
    ui.extra{end+1} = ui.type;
    hymos_visual_style(ui.type, 'FontSize', 18, 'FontWeight', 'bold', ...
                       'FontColor', [0.15 0.31 0.58]);
    for i = 1:8
        y = 586-(i-1)*62;
        ui.labels{i} = uilabel(fig, 'Text', '', 'Position', [280 y+7 190 24]);
        ui.extra{end+1} = ui.labels{i};
        ui.edits{i} = uieditfield(fig, 'text', 'Value', '', 'Position', [485 y 210 38]);
        ui.extra{end+1} = ui.edits{i};
        ui.choices{i} = uidropdown(fig, 'Items', {'single', 'nmg'}, ...
                                  'Value', 'single', 'Position', [485 y 210 38]);
        ui.extra{end+1} = ui.choices{i};
        ui.hints{i} = uilabel(fig, 'Text', '', 'Position', [715 y+7 405 24]);
        ui.extra{end+1} = ui.hints{i};
        hymos_visual_style(ui.labels{i}, 'FontSize', 18, 'FontColor', [0.15 0.21 0.32]);
        hymos_visual_style(ui.edits{i}, 'FontSize', 19);
        hymos_visual_style(ui.choices{i}, 'FontSize', 18);
        hymos_visual_style(ui.hints{i}, 'FontSize', 14, 'FontColor', [0.43 0.49 0.59]);
    end
    % Bind callbacks only after all controls exist.  The host captures
    % callback arguments by value, so binding earlier leaves later cells empty.
    for i = 1:8
        set(ui.edits{i}, 'ValueChangedFcn', {@hymos_form_value_changed, fig, ui});
        set(ui.choices{i}, 'ValueChangedFcn', {@hymos_form_value_changed, fig, ui});
    end
catch err
    ui.form_available = false;
    disp('分类表单控件不可用，已回退到完整文本编辑。');
    disp(err);
end
% Optional guidance must never prevent the main form from being created on
% hosts that implement only a subset of modern UI control properties.
try
    ui.advanced_help = uitextarea(fig, 'Position', [260 555 875 87], 'Value', { ...
        '优先使用左侧表单；这里只在需要查看或修改完整配置时使用。', ...
        '格式：[System] 是分组；Key = Value 是设置；# 开头是说明。', ...
        '修改后先点“检查参数”，通过后再“应用并关闭”。'});
    set(ui.advanced_help, 'Enable', 'off');
    hymos_visual_style(ui.advanced_help, 'FontSize', 14, ...
        'BackgroundColor', [0.94 0.97 1], 'FontColor', [0.20 0.29 0.42]);
catch
    ui.advanced_help = [];
end
% UserData is plain draft data, never a Couette task or task-containing object.
data.text = state.text; data.path = state.path; data.default_path = state.default_path;
data.page = 5; data.rows = []; data.values = {}; data.displayed_text = '';
set(fig, 'UserData', data);
for i = 1:5
    set(ui.nav{i}, 'ButtonPushedFcn', {@hymos_form_navigate, i, fig, ui});
    if ~ui.form_available && i < 5
        set(ui.nav{i}, 'Enable', 'off');
    end
end
set(ui.editor, 'ValueChangedFcn', {@hymos_form_dirty, ui.message});
set(load_button, 'ButtonPushedFcn', {@hymos_form_action, 'load', fig, ui, template_path});
set(save_button, 'ButtonPushedFcn', {@hymos_form_action, 'save', fig, ui, template_path});
set(save_as_button, 'ButtonPushedFcn', {@hymos_form_action, 'save_as', fig, ui, template_path});
set(default_button, 'ButtonPushedFcn', {@hymos_form_action, 'default', fig, ui, template_path});
set(validate_button, 'ButtonPushedFcn', {@hymos_form_action, 'validate', fig, ui, template_path});
set(apply_button, 'ButtonPushedFcn', {@hymos_form_action, 'apply', fig, ui, template_path});
try
    if ui.form_available
        hymos_form_render(fig, ui, 1);
    else
        hymos_form_render(fig, ui, 5);
    end
catch err
    hymos_form_render(fig, ui, 5);
    disp(err);
end
set(ui.message, 'Text', [state.message ' 当前参数待校验。']);
end

function lines = hymos_form_lines(text)
lines = strsplit(regexprep(char(text), '\r\n?', '\n'), newline);
end

function text = hymos_form_read_editor(editor)
value = get(editor, 'Value');
if iscell(value)
    text = strjoin(value, newline);
else
    text = char(value);
end
end

function hymos_form_dirty(src, event, message)
set(message, 'Text', '参数已修改，待校验；尚未保存或应用。');
end

function hymos_form_value_changed(src, event, fig, ui)
data = get(fig, 'UserData');
for i = 1:length(data.rows)
    if isequal(src, ui.edits{i}) || isequal(src, ui.choices{i})
        key = ui.schema.rows{data.rows(i), 2};
        if data.page == 1 && hymos_form_is_wall_key(key)
            value = char(get(ui.edits{i}, 'Value'));
            data.text = char(hymos_ui_config(data.text, key, value));
            set(fig, 'UserData', data);
            set(ui.type, 'Text', ['当前工况：' ...
                hymos_form_type_text(hymos_form_case_label(data.text))]);
        end
        if strcmp(key, 'Solver/SolverType')
            hymos_form_set_multigrid_nav(ui, char(get(ui.choices{i}, 'Value')));
        end
        break;
    end
end
hymos_form_dirty(src, event, ui.message);
end

function tf = hymos_form_is_wall_key(key)
tf = strcmp(key, 'LeftWall/Temperature') || strcmp(key, 'RightWall/Temperature') || ...
     strcmp(key, 'LeftWall/TangentialVelocity') || strcmp(key, 'RightWall/TangentialVelocity');
end

function label = hymos_form_case_label(text)
left_t = str2double(char(hymos_ui_config(text, 'LeftWall/Temperature')));
right_t = str2double(char(hymos_ui_config(text, 'RightWall/Temperature')));
left_u = str2double(char(hymos_ui_config(text, 'LeftWall/TangentialVelocity')));
right_u = str2double(char(hymos_ui_config(text, 'RightWall/TangentialVelocity')));
if ~isfinite(left_t) || ~isfinite(right_t) || ~isfinite(left_u) || ~isfinite(right_u)
    label = 'Custom';
    return;
end
same_temperature = hymos_form_same_value(left_t, right_t);
same_velocity = hymos_form_same_value(left_u, right_u);
if same_temperature && ~same_velocity
    label = 'Couette';
elseif ~same_temperature && same_velocity
    label = 'Fourier';
elseif ~same_temperature && ~same_velocity
    label = 'Couette–Fourier';
else
    label = 'Equilibrium';
end
end

function text = hymos_form_type_text(label)
if strcmp(label, 'Custom')
    text = '参数无效';
elseif strcmp(label, 'Equilibrium')
    text = '平衡态 / 无壁面驱动';
else
    text = label;
end
end

function tf = hymos_form_same_value(left, right)
tf = abs(left - right) <= 1e-12 * max([1, abs(left), abs(right)]);
end

function data = hymos_form_collect(fig, ui)
data = get(fig, 'UserData');
keys = {}; values = {};
if data.page == 5
    values = {hymos_form_read_editor(ui.editor)};
else
    for i = 1:length(data.rows)
        row = data.rows(i);
        keys{end+1} = ui.schema.rows{row, 2};
        if isempty(ui.schema.rows{row, 5})
            value = char(get(ui.edits{i}, 'Value'));
        else
            value = char(get(ui.choices{i}, 'Value'));
            if strcmp(value, '(未设置)')
                value = '';
            end
        end
        values{end+1} = value;
    end
end
data = hymos_form_sync(data, keys, values);
set(fig, 'UserData', data);
end

function hymos_form_render(fig, ui, page)
data = get(fig, 'UserData');
rows = []; values = {};
if page ~= 5
    for row = 1:size(ui.schema.rows, 1)
        if ui.schema.rows{row, 1} == page
            rows(end+1) = row;
            value = char(hymos_ui_config(data.text, ui.schema.rows{row, 2}));
            if strcmp(ui.schema.rows{row, 2}, 'MultiGridMethod/CycleType') && isempty(value)
                value = char(hymos_ui_config(data.text, 'MultiGridMethod/VW_CYCLE'));
            end
            values{end+1} = value;
        end
    end
end
% Parse before hiding the current page. Invalid advanced text stays editable.
for i = 1:length(ui.extra)
    set(ui.extra{i}, 'Visible', 'off');
end
if ~isempty(ui.advanced_help)
    set(ui.advanced_help, 'Visible', 'off');
end
set(ui.editor, 'Visible', 'off');
set(ui.title, 'Text', '一维平行板 BGK 算例');
set(ui.section, 'Text', ['   ' ui.schema.groups{page}]);
set(ui.note, 'Text', '只修改当前分类；其余配置及注释保留在完整文本中。');
for i = 1:5
    label = ui.schema.groups{i};
    if i == page
        hymos_visual_style(ui.nav{i}, 'BackgroundColor', [0.83 0.91 1], ...
            'FontColor', [0.06 0.31 0.83], 'FontWeight', 'bold');
    else
        hymos_visual_style(ui.nav{i}, 'BackgroundColor', [0.94 0.96 0.99], ...
            'FontColor', [0.28 0.34 0.44], 'FontWeight', 'normal');
    end
    set(ui.nav{i}, 'Text', label);
end
hymos_form_set_multigrid_nav(ui, char(hymos_ui_config(data.text, 'Solver/SolverType')));
if page == 5
    editor_position = [260 180 875 462];
    if ~isempty(ui.advanced_help)
        set(ui.advanced_help, 'Visible', 'on');
        editor_position = [260 180 875 365];
    end
    set(ui.editor, 'Position', editor_position, ...
        'Value', hymos_form_lines(data.text), 'Visible', 'on');
    data.displayed_text = hymos_form_read_editor(ui.editor);
    set(ui.note, 'Text', '完整 TXT/INI；修改后请先检查参数。');
else
    enabled = 'on';
    if page == 2
        mesh_type = char(hymos_ui_config(data.text, 'Mesh/Type'));
        set(ui.note, 'Text', ['一维均匀网格；仅支持 Mesh/Type=0。当前 Type=' mesh_type]);
    elseif page == 4
        solver = char(hymos_ui_config(data.text, 'Solver/SolverType'));
        if strcmp(solver, 'single')
            enabled = 'off';
            set(ui.note, 'Text', 'single 不使用多重网格。');
        else
            set(ui.note, 'Text', '多重网格参数。');
        end
    end
    if page == 1
        set(ui.type, 'Text', ['当前工况：' hymos_form_type_text(hymos_form_case_label(data.text))], ...
            'Visible', 'on');
        set(ui.note, 'Text', '当前工况由左右壁温和切向壁速自动识别。');
    end
    % Place controls evenly through the page while keeping the case label clear.
    spacing = 62; first_y = 578;
    if page == 1 && length(rows) >= 5
        spacing = 62; first_y = 530;
    elseif length(rows) >= 7
        spacing = 68; first_y = 590;
    elseif length(rows) == 6
        spacing = 72; first_y = 560;
    elseif length(rows) <= 4
        spacing = 112; first_y = 550;
    elseif length(rows) == 5
        spacing = 86; first_y = 560;
    end
    for i = 1:length(rows)
        y = first_y-(i-1)*spacing;
        set(ui.labels{i}, 'Position', [280 y+8 195 24]);
        set(ui.edits{i}, 'Position', [485 y 210 40]);
        set(ui.choices{i}, 'Position', [485 y 210 40]);
        set(ui.hints{i}, 'Position', [715 y+8 405 24]);
        row = rows(i);
        set(ui.labels{i}, 'Text', ui.schema.rows{row, 3}, 'Visible', 'on');
        set(ui.hints{i}, 'Text', ui.schema.rows{row, 4}, 'Visible', 'on');
        items = ui.schema.rows{row, 5};
        if isempty(items)
            set(ui.edits{i}, 'Value', values{i}, 'Enable', enabled, 'Visible', 'on');
        else
            selected = values{i};
            if isempty(selected)
                selected = '(未设置)';
            end
            present = false;
            for k = 1:length(items)
                if strcmp(items{k}, selected)
                    present = true;
                end
            end
            if ~present
                items{end+1} = selected;
            end
            set(ui.choices{i}, 'Items', items);
            set(ui.choices{i}, 'Value', selected, 'Enable', enabled, 'Visible', 'on');
        end
    end
end
data.page = page; data.rows = rows; data.values = values;
set(fig, 'UserData', data);
set(ui.path, 'Text', ['当前文件：' data.path]);
set(ui.message, 'Text', '参数待校验；保存与应用均使用原配置校验。');
end

function hymos_form_set_multigrid_nav(ui, solver)
if strcmpi(strtrim(solver), 'single')
    set(ui.nav{4}, 'Enable', 'off');
    hymos_visual_style(ui.nav{4}, 'BackgroundColor', [0.94 0.96 0.99], ...
        'FontColor', [0.62 0.66 0.72], 'FontWeight', 'normal');
else
    set(ui.nav{4}, 'Enable', 'on');
end
end

function hymos_form_navigate(src, event, page, fig, ui)
try
    data = hymos_form_collect(fig, ui);
catch err
    set(ui.message, 'Text', '当前输入不能转为配置，请修正后再切换。');
    disp(err);
    return;
end
try
    hymos_form_render(fig, ui, page);
catch err
    % Covers both malformed text and partially unsupported form controls.
    hymos_form_render(fig, ui, 5);
    set(ui.message, 'Text', '不能显示该表单；原文保留在完整文本中。详情见命令窗口。');
    disp(err);
end
end

function hymos_form_action(src, event, action, fig, ui, template_path)
try
    if strcmp(action, 'load')
        [file, location] = uigetfile({'*.txt;*.ini', 'HyMoS 参数文件'}, '加载参数');
        if isequal(file, 0)
            return;
        end
        text = fileread(fullfile(location, file));
        data = get(fig, 'UserData');
        data.text = text;
        try
            data.path = char(hymos_ui_file(fullfile(location, file)));
        catch
            data.path = data.default_path;
            disp('源文件不可写回；保存目标改为当前工作目录。');
        end
        set(fig, 'UserData', data);
        hymos_form_render(fig, ui, 5);
        set(ui.message, 'Text', '已加载完整文本，尚未应用；可从左侧切换分类。');
        return;
    elseif strcmp(action, 'default')
        text = fileread(char(template_path));
        check = hymos_validate_text('couette', text);
        if ~check.valid
            error(char(check.message));
        end
        data = get(fig, 'UserData');
        data.text = text; data.path = data.default_path;
        set(fig, 'UserData', data);
        hymos_form_render(fig, ui, 5);
        set(ui.message, 'Text', '已恢复默认草稿；尚未保存或应用，运行任务不受影响。');
        return;
    end
    data = hymos_form_collect(fig, ui);
    if strcmp(action, 'validate')
        check = hymos_validate_text('couette', data.text);
        if check.valid
            set(ui.message, 'Text', '原配置输入校验通过（不代表所有组合收敛）；尚未应用。');
        else
            set(ui.message, 'Text', '配置无效；原文保留，详情见命令窗口。');
            disp(char(check.message));
        end
    elseif strcmp(action, 'save') || strcmp(action, 'save_as')
        path = data.path;
        if strcmp(action, 'save_as')
            [file, location] = uiputfile('*.txt', '参数另存为', path);
            if isequal(file, 0)
                return;
            end
            path = fullfile(location, file);
        end
        data.path = char(hymos_ui_file(path, data.text));
        set(fig, 'UserData', data);
        set(ui.path, 'Text', ['当前文件：' data.path]);
        set(ui.message, 'Text', '参数文件已保存；会话参数尚未应用。');
        disp(['HyMoS 参数已保存：' data.path]);
    elseif strcmp(action, 'apply')
        path = hymos_editor_apply(data.text, data.path, get(ui.save_on_apply, 'Value'));
        disp('HyMoS 参数已应用。运行并监视：hymos_run();');
        try
            close(fig);
        catch err
            set(ui.message, 'Text', '参数已应用，但窗口未关闭；请手动关闭。');
            disp(err);
        end
    end
catch err
    set(ui.message, 'Text', '操作未完成，窗口和草稿保留。详情见命令窗口。');
    disp(err);
end
end

function hymos_visual_style(control, varargin)
% Cosmetic properties must never disable parameter editing on older hosts.
% Apply separately so an unsupported color does not suppress a supported font.
persistent reported;
for k = 1:2:length(varargin)
    try
        set(control, varargin{k}, varargin{k+1});
    catch
        if isempty(reported)
            disp('部分外观属性不受宿主支持，保留默认样式；参数功能不受影响。');
            reported = true;
        end
    end
end
end
