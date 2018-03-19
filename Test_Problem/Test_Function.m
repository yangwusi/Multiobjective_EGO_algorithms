function [design_space,ref_point] = Test_Function(name, num_obj, num_vari)
%--------------------------------------------------------------------------
% ZDT with 2 objectives
switch name
    case {'ZDT1', 'ZDT2', 'ZDT3', 'ZDT5', 'ZDT6'}
        design_space=[zeros(1,num_vari);ones(1,num_vari)];
        ref_point=11*ones(1,2);
    case 'ZDT4'
        design_space=[0,-5*ones(1,num_vari-1);1,5*ones(1,num_vari-1)];
        ref_point=11*ones(1,2);
    otherwise % the DTLZ test problems
        design_space=[zeros(1,num_vari);ones(1,num_vari)];
        if num_obj==2
            ref_point=10*ones(1,num_obj);
        else
            if str2double(name(5))==2 || str2double(name(5))==5
                ref_point=2.5*ones(1,num_obj);
            elseif str2double(name(5))==7
                if num_obj==3
                    ref_point=30*ones(1,num_obj);
                elseif num_obj==4
                    ref_point=50*ones(1,num_obj);
                elseif num_obj==5
                    ref_point=60*ones(1,num_obj);
                elseif num_obj==6
                    ref_point=70*ones(1,num_obj);
                end
            end
        end
end
end
