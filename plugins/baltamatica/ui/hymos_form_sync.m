function data = hymos_form_sync(data, keys, values)
% Transactional draft update shared by GUI actions and host regression tests.
% Advance the displayed baseline after every collection, including Save and
% Validate; otherwise changing a field back to its earlier value loses edits.
if data.page == 5
    visible = values{1};
    if ~strcmp(visible, data.displayed_text)
        data.text = visible;
    end
    data.displayed_text = visible;
else
    candidate = data.text;
    for i = 1:length(keys)
        if ~strcmp(values{i}, data.values{i})
            candidate = char(hymos_ui_config(candidate, keys{i}, values{i}));
        end
    end
    data.text = candidate;
    data.values = values;
end
end
