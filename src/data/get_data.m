function data = get_data(type, model, modes, varargin)
    
switch lower(type)
    case 'foam_leap'
        data = struct();
        data.type  = "foam_leap";
        data.modes = upper(string(modes));
        for i = 1:numel(modes)
            mode = upper(modes{i});
            data.(mode) = load_foam_leap(mode);
        end
    case 'foam_turbo'
        data = struct();
        data.type  = "foam_turbo";
        data.modes = upper(string(modes));
        for i = 1:numel(modes)
            mode = upper(modes{i});
            data.(mode) = load_foam_turbo(mode);
        end

    otherwise
        error('Unknown data type');
end
end
