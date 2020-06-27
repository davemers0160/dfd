function data = read_dfd_images(input)
    %custom inputfile for reading in dfd infous and out of focus images.
    % this assumes that the input is the infocus image and then the our of
    % focus image needs to be read in too.
    
    df = '_lin_0.32_2.88.png';
    im_f = imread(input);
    
    file_base = strsplit(input,'.png');
    im_b = imread(strcat(file_base{1},df));
    data = cat(3,im_f,im_b);
    
end