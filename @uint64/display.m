% Overloaded display function for uint64 
% (elements of GF(q^22)
%
% Returns signed value

function display(X)
if isequal(get(0,'FormatSpacing'),'compact')
   disp([inputname(1) ' =']);
   disp(convert(X))
else
   disp(' ')
   disp([inputname(1) ' =']);
   disp(' ');
   disp(convert(X))
end