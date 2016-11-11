function stab_sites = load_stab_sites(file)
    fileID = fopen(file,'r');

    stab_sites = textscan(fileID,'%s');
    stab_sites = stab_sites{1};
    
    fclose(fileID);
end