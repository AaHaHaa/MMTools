function A = catstruct(varargin)
% CATSTRUCT   Concatenate or merge structures with different fieldnames
%   X = CATSTRUCT(S1,S2,S3,...) merges the structures S1, S2, S3 ...
%   into one new structure X. X contains all fields present in the various
%   structures. An example:
%
%     A.name = 'Me' ;
%     B.income = 99999 ;
%     X = catstruct(A,B) 
%     % -> X.name = 'Me' ;
%     %    X.income = 99999 ;
%
%   If a fieldname is not unique among structures (i.e., a fieldname is
%   present in more than one structure), only the value from the last
%   structure with this field is used. In this case, the fields are 
%   alphabetically sorted. A warning is issued as well. An axample:
%
%     S1.name = 'Me' ;
%     S2.age  = 20 ; S3.age  = 30 ; S4.age  = 40 ;
%     S5.honest = false ;
%     Y = catstruct(S1,S2,S3,S4,S5) % use value from S4
%
%   The inputs can be array of structures. All structures should have the
%   same size. An example:
%
%     C(1).bb = 1 ; C(2).bb = 2 ;
%     D(1).aa = 3 ; D(2).aa = 4 ;
%     CD = catstruct(C,D) % CD is a 1x2 structure array with fields bb and aa
%
%   The last input can be the string 'sorted'. In this case,
%   CATSTRUCT(S1,S2, ..., 'sorted') will sort the fieldnames alphabetically. 
%   To sort the fieldnames of a structure A, you could use
%   CATSTRUCT(A,'sorted') but I recommend ORDERFIELDS for doing that.
%
%   When there is nothing to concatenate, the result will be an empty
%   struct (0x0 struct array with no fields).
%
%   NOTE: To concatenate similar arrays of structs, you can use simple
%   concatenation: 
%     A = dir('*.mat') ; B = dir('*.m') ; C = [A ; B] ;
%
%   NOTE: This function relies on unique. Matlab changed the behavior of
%   its set functions since 2013a, so this might cause some backward
%   compatibility issues when dulpicated fieldnames are found.
%
%   See also CAT, STRUCT, FIELDNAMES, STRUCT2CELL, ORDERFIELDS
%
% version 4.1 (feb 2015), tested in R2014a
% (c) Jos van der Geest
% email: jos@jasen.nl
%
% =========================================================================
% ----------------------- Additional functionality ------------------------
%   Deal with the multiple-layer structures, a structure inside another
%   structure.
% =========================================================================
%   The last one or two input arguments should be a structure,
%   "sorted" or "overwrite".
%   The strings don't have to be in a certain order.
%
%       E.g.
%           s = catstruct(s1,s2,'sorted');
%           s = catstruct(s1,s2,'overwrite');
%           s = catstruct(s1,s2,'overwrite','sorted');
%           s = catstruct(s1,s2,'sorted','overwrite');
%
%   "sorted" is the same as the above.
%   "overwrite":
%       (1) If there's no "overwrite", it'll run "catstruct" again for the
%       structures under a structure.
%       (2) If "overwrite", it works as the original "catstruct", that is,
%           the structure under a structure might be overwritten if there's
%           the same field.
%
%       For example,
%
%           s1.a = 1; s1.b.ba = 2; s1.b.bb = 3; s1.c = 4; % s1.b is a structure
%           s2.a = 5; s2.b.ba = 6; s2.b.bc = 7; s2.d = 8; % s2.b is a structure
%
%           s = catstruct(s1,s2,'overwrite'); % original catstruct
%
%           Then, s: a = 5;
%                    b.ba = 6;    b.bc = 7; % s2.b overwrites s1.b
%                    c = 4;
%                    d = 8;
%
%           s = catstruct(s1,s2);
%           
%           Then, s: a = 5;
%                    b.ba = 6;    b.bb = 3;    b.bc = 7; % catstruct(s1.b,s2.b)
%                    c = 4;
%                    d = 8;
%
% version 5.0 (jan 2018), tested in R2017a, R2017b
% (c) Yi-Hao Chen
% email: yc2368@cornell.edu
%
% History
% Created in 2005
% Revisions
%   2.0 (sep 2007) removed bug when dealing with fields containing cell
%                  arrays (Thanks to Rene Willemink)
%   2.1 (sep 2008) added warning and error identifiers
%   2.2 (oct 2008) fixed error when dealing with empty structs (thanks to
%                  Lars Barring)
%   3.0 (mar 2013) fixed problem when the inputs were array of structures
%                  (thanks to Tor Inge Birkenes).
%                  Rephrased the help section as well.
%   4.0 (dec 2013) fixed problem with unique due to version differences in
%                  ML. Unique(...,'last') is no longer the deafult.
%                  (thanks to Isabel P)
%   4.1 (feb 2015) fixed warning with narginchk
%   5.0 (jan 2018) added functionality of merging multiple layers of
%                  structures

narginchk(1,Inf);
N = nargin;

sorted = false;
overwrite = false;
while ~isstruct(varargin{N})
    switch varargin{N}
        case 'sorted'
            if sorted
                error('catstruct:InvalidArgument','"sorted" is repeated in the arguments.');
            end
            sorted = true;
        case 'overwrite'
            if overwrite
                error('catstruct:InvalidArgument','"overwrite" is repeated in the arguments.');
            end
            overwrite = true;
        otherwise
            error('catstruct:InvalidArgument','Last argument should be a structure, or the string "sorted", "overwrite".') ;
    end
    N = N-1;
    narginchk(nargin-N+1,Inf);
end

sz0 = []; % used to check that all inputs have the same size

% used to check for a few trivial cases
NonEmptyInputs = false(N,1);
NonEmptyInputsN = 0;

% used to collect the fieldnames and the inputs
FN = cell(N,1);
VAL = cell(N,1);

% parse the inputs
for ii=1:N
    X = varargin{ii};
    if ~isstruct(X)
        error('catstruct:InvalidArgument',['Argument #' num2str(ii) ' is not a structure.']);
    end
    
    if ~isempty(X)
        % empty structs are ignored
        if ii > 1 && ~isempty(sz0)
            if ~isequal(size(X), sz0)
                error('catstruct:UnequalSizes','All structures should have the same size.');
            end
        else
            sz0 = size(X) ;
        end
        NonEmptyInputsN = NonEmptyInputsN + 1;
        NonEmptyInputs(ii) = true;
        FN{ii} = fieldnames(X);
        VAL{ii} = struct2cell(X);
    end
end

switch NonEmptyInputsN
    case 0
        % all structures were empty
        A = struct([]) ;
    case 1
        % there was only one non-empty structure
        A = varargin{NonEmptyInputs};
        if sorted
            A = orderfields(A);
        end
    otherwise
        % there is actually something to concatenate
        FN = cat(1,FN{:});    
        VAL = cat(1,VAL{:});    
        FN = squeeze(FN);
        VAL = squeeze(VAL);

        FN_origin = FN;

        [UFN,ind] = unique(FN, 'last');
        % If this line errors, due to your matlab version not having UNIQUE
        % accept the 'last' input, use the following line instead
        % [UFN,ind] = unique(FN) ; % earlier ML versions, like 6.5

        if numel(UFN) ~= numel(FN)
            warning('catstruct:DuplicatesFound','Fieldnames are not unique between structures.');
            sorted = true;
        end

        if sorted
            VAL = VAL(ind,:);
            FN = FN(ind,:);
        end

        A = cell2struct(VAL, FN, 1);
        A = reshape(A, sz0) ; % reshape into original format

        % Added by Yi-Hao Chen:
        % Deal with the multiple-layer structures, a structure under another structure.
        if ~overwrite
            [UFN,~, FN_indices] = unique(FN_origin,'stable');
            [n, bin] = histc(FN_indices, unique(FN_indices));
            multiple = find(n > 1);
            index = find(ismember(bin, multiple));
            sameField = unique(FN_indices(index));
            if ~isempty(sameField)
                for sameField_in_loop = sameField'
                    for structarrayIndex = 1:length(X)
                        k = 1;
                        struct_under_struct = cell(length(find(FN_indices(index)==sameField_in_loop)),1);
                        for ii = 1:N
                            X = varargin{ii};
                            FN_each = fieldnames(X);
                            if any(strcmp(FN_each,UFN{sameField_in_loop}))
                                if ~isstruct(X(structarrayIndex).(UFN{sameField_in_loop}))
                                    break;
                                else
                                    struct_under_struct{k} = X(structarrayIndex).(UFN{sameField_in_loop});
                                    k = k + 1;
                                end
                            end
                        end
                        if ii ~= 1
                            if sorted
                                merged_struct = catstruct(struct_under_struct{:},'sorted');
                            else
                                merged_struct = catstruct(struct_under_struct{:});
                            end
                            A(structarrayIndex).(UFN{sameField_in_loop}) = merged_struct;
                        end
                    end
                end
            end
        end
    
end



