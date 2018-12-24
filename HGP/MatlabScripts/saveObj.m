function [  ] = saveObj( path, name, mesh )

obj.vertices = mesh.V;
obj.vertices_texture = mesh.UV;
objects.data.vertices = mesh.F;
objects.data.texture = mesh.halfEdges;
obj.objects = objects;
obj.objects.type='f';

write_wobj(obj,[path '\\' name])

end

