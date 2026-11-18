function Res = mergestruct(A,B)
% MERGESTRUCT  Recursively merge the fields of structure B into structure A
%
% Fields of B overwrite those of A; struct-valued fields present in both are
% merged recursively. An empty or yaml.Null field in A is replaced by B.
%
% Use as:
%   Res = mergestruct(A, B)
%
% Example:
%   A.field1=1;
%   A.field2.subfield1=1;
%   A.field2.subfield2=2;
%
%   B.field1=1;
%   B.field2.subfield1=10;
%   B.field2.subfield3=30;
%   B.field3.subfield1=1;
%
%   C=mergestruct(A,B);
%

Res=[];
if nargin>0
    Res=A;
end;
if nargin==1 || isstruct(B)==0
    return;
end;
  fnb=fieldnames(B);

  for i=1:length(fnb)
     s=char(fnb(i));
     oldfield=[];
     if (isfield(A,s))
         oldfield=getfield(A,s);
     end
     newfield=getfield(B,s);
     % treat yaml.Null as empty so YAML-loaded null values are always overwritten
     if isempty(oldfield) || isstruct(newfield)==0 || isa(oldfield, 'yaml.Null')
       Res=setfield(Res,s,newfield);
     else
       Res=setfield(Res,s,mergestruct(oldfield, newfield));
     end
  end

end
