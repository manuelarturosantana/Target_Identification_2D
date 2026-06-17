function model = attach_glrt_classification_map(model, eval_data, library_data, true_to_library_map)
%ATTACH_GLRT_CLASSIFICATION_MAP  Attach optional truth metadata to GLRT model.
%
%   model = ATTACH_GLRT_CLASSIFICATION_MAP(model, eval_data, library_data,
%   true_to_library_map) attaches metadata used by classify_glrt_data when
%   the classified objects and pole-library objects are not the same list.
%
%   true_to_library_map(i,j) is true if evaluation object i should be
%   counted correct when classified as library/model object j.

    neval = numel(eval_data);
    nlib = numel(library_data);

    if numel(model) ~= nlib
        error('attach_glrt_classification_map:ModelLibraryMismatch', ...
              'numel(model) must equal numel(library_data).');
    end

    if ~isequal(size(true_to_library_map), [neval, nlib])
        error('attach_glrt_classification_map:BadMapSize', ...
              'true_to_library_map must have size %d-by-%d.', neval, nlib);
    end

    true_to_library_map = logical(true_to_library_map);

    if any(sum(true_to_library_map, 2) == 0)
        bad = find(sum(true_to_library_map, 2) == 0, 1, 'first');
        error('attach_glrt_classification_map:EmptyTruthRow', ...
              'Evaluation object %d has no acceptable library prediction.', bad);
    end

    cls = struct();
    cls.true_to_library_map = true_to_library_map;
    cls.library_data = library_data;
    cls.library_names = strings(nlib, 1);
    cls.eval_names = strings(neval, 1);

    for j = 1:nlib
        cls.library_names(j) = library_data(j).name;
    end

    for i = 1:neval
        cls.eval_names(i) = eval_data(i).name;
    end

    model(1).classification = cls;
end