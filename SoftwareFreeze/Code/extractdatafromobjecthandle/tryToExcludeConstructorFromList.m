function out = tryToExcludeConstructorFromList(in, obj)
    classNameWithoutFullPackageAddress = regexp(class(obj), '[^.]*$', 'match');
    out = in(~strcmp(in, classNameWithoutFullPackageAddress));
end