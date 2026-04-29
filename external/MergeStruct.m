function Res = mergestruct(A,B)
%% Recursively merges fields and subfields of structures A and B to result structure Res
% Simple recursive algorithm merges fields and subfields of two structures
%   Example:
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
     if isempty(oldfield) || isstruct(newfield)==0
       Res=setfield(Res,s,newfield);     
     else
       Res=setfield(Res,s,MergeStruct(oldfield, newfield));  
     end    
  end    

end

