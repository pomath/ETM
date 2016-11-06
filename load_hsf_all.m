function A = load_hsf_all(epoch_i,epoch_j,stnm)
    
    A = zeros(1,size(stnm,1));
    files = dir('stations/H_*');
    
    for i = 1:size(files,1)
        pos = strmatch(files(i).name(3:6),stnm);
        if pos ~= 0
            Ha = load(['stations/' files(i).name]);

            for j = 1:size(Ha,1)
                if Ha(j,1) >= epoch_i & Ha(j,1) <= epoch_j
                    A(pos) = 1;
                end
            end
        end
    end
end