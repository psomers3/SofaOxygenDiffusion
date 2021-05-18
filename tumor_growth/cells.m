classdef cells  
    properties
        altePos
        pos
        radius
        teilbar=false
        nachbarn=cells.empty
        F_rep=0
        F_adh=0
        f=0
        dist
        F=[0 0 0]
        v=[0 0 0]
        fibreForce=[0 0 0]
        maxradius=5e-6
        phenotype=0; %0 normoxic, 1 hypoxic, 2 necrotic
        oxygen=60
        ID
    end
    methods
        function obj = grow(obj)
            %growth of a cell
            if obj.radius<obj.maxradius
                obj.radius=obj.radius+0.1e-6;
            end
            if obj.radius> obj.maxradius
                obj.radius= obj.maxradius;
            end
            if obj.radius>0.99* obj.maxradius
        		obj.teilbar=true;
            end
        end
        function obj = move3D(obj)
            %movement of a cell
            if obj.phenotype==0
                ar=2*10^-3*10^-6;
            end
            if obj.phenotype==1
                ar=2*10^-3*10^-6*10;
            end
            if obj.phenotype==2
                ar=0;
            end
            timestep=1;
            gamma=24;
            obj.pos=obj.pos+timestep/gamma*(ar*[randn randn randn]+obj.F);
        end
    end
end




    
    