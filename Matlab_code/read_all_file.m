function [ allnames, len ] = read_all_file( path,current_path )
    cd(path);
    list=dir('*.mat');
    name={list.name}; 
    allnames=natsort(name);
    [k,len]=size(allnames); 
    cd(current_path);
end




 
 