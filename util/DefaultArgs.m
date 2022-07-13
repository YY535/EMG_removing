function varargout = DefaultArgs(Ins, InsC)
varargout = cell(length(InsC));
lIns = length(Ins);
for k = 1:lIns
    if ~isempty(Ins{k})
        varargout{k} = Ins{k};
    else
        varargout{k} = InsC{k};
    end
        
end
for k = (lIns+1):length(InsC)
    varargout{k} = InsC{k};
end
