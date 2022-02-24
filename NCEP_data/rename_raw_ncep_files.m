ccc


years = [2002,2003];

d = dir;
for j = 1:length(years)
    for i = 1:length(d)
        fn = d(i).name;
        if length(fn)>4 && contains(fn,string(years(j)))
            dot_vec = strfind(fn,'.');
            fn_dash = fn;
            for k = 1:length(dot_vec)-1
                fn_dash(dot_vec(k)) = '_';
            end
            copyfile(fn,fn_dash)
        end
    end
end
    
    
    
    
    
    
    
    
    
    
    
    
    
