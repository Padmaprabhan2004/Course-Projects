%A.PADMAPRABHAN-ME22BTECCH11001
%QUESTION 5 - Part (d) Discretization and global assembly of beam bending
%problem. Assumed 2,4 and 8 elements and assembled global stiffness and
%force vector for ech case.
function BeamBendingAssembly()    
    L = 10; 
    EI = 3;  
    q = 1;   
    no_of_ele_arr = [2, 4, 8];
    for no_of_elements = no_of_ele_arr
        le = L / no_of_elements;  
        K_glob = zeros(2*no_of_elements + 2); 
        F_glob = zeros(2*no_of_elements + 2, 1); 
        K_local = (EI / le^3) * [[12,6*le,-12,6*le];
            [6*le,4*le^2,-6*le,2*le^2];
            [-12,-6*le,12,-6*le];
            [6*le,2*le^2,-6*le,4*le^2]];
        for i = 1:no_of_elements
            glob_indices = 2*(i-1) + (1:4);
            K_glob(glob_indices, glob_indices) = K_glob(glob_indices, glob_indices) + K_local;
            F_glob(glob_indices) = F_glob(glob_indices) + le/2 * [q; q*le/6; q; -q*le/6];
        end
       K_glob(1,:)=0;
       K_glob(:,1)=0;
       K_glob(1,1)=1;
       F_glob(1)=0;
       K_glob(end,:)=0;
       K_glob(:,end)=0;
       K_glob(end,end)=1;
       F_glob(end)=0;
       %global assembly
       w = K_glob \ F_glob;
       fprintf('Deflections for %d elements:\n', no_of_elements);
       disp(w');
    end
end