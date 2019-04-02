function onew = choseSubset(obj, subset)

%% Function that selects subset from object in the ZBraingrid class.
%
%  Takes in input the properties of the object, and find comments from each
%  datasets that include the subset string. This is like doing an indexing,
%  but instead of giving indexes, an expression is used to discriminate
%  subsets.
% 
%
%% Inputs:
%
%  --obj: references the object this methods is attached to.
%  --subset: way of selecting the subset to be plot. String, which will
%    trigger a research in the comment property of the ZBraingrid object.
%
%
%% Output:
%
%  --onew: new ZBraingrid object, sub-object from obj.



    %% Converting string research to array of subsets:

    keepsub = [];
    for i = 1:length(obj.comments)
        if ~isempty(strfind(obj.comments{i}, subset))
            keepsub = [keepsub, i];
        end
    end
    if isempty(keepsub)
        error('Subset to plot is empty.')
    end
    
    
    
    %% Using indexing method to obtain onew:
    
    Sin = struct;
    Sin.type = '()';
    Sin.subs = {keepsub};
    onew = subsref(obj, Sin);
    
    
    
end