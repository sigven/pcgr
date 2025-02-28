function [sum ] = hashseq(seq) 
    sum = 0;
    len = length(seq);
    four = [1;4;16;64;256;1024;4096];
    for i = 1:len
        sum = sum + str2num(seq(i)) * four(len - i +1) ;
    end
    