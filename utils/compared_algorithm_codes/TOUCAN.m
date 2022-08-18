function Xhat=TOUCAN(MissM,Mask,I,option)

if ndims(MissM)~=3
    error('Ndims not 3');
end
dim = size(MissM);

MissMt = permute(MissM,[3 2 1]);
Maskt = permute(Mask,[3 2 1]);

obj = py.importlib.import_module('toucan_mod');
py.importlib.reload(obj);
output=py.toucan_mod.toucan(MissMt(:),Maskt(:),uint16(dim));
data = double(py.array.array('d',py.numpy.nditer(output{1})));
Xhat=permute(reshape(data,dim(2),dim(3),dim(1)),[3 1 2]);

end