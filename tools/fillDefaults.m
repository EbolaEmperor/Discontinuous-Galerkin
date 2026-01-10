function opts = fillDefaults(opts, defaults)
    fields = fieldnames(defaults);
    for k = 1:numel(fields)
        f = fields{k};
        if ~isfield(opts, f) || isempty(opts.(f))
            opts.(f) = defaults.(f);
        end
    end
end