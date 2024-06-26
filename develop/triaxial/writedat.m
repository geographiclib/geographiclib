function writedat(dat,file)
  fid = fopen(file, 'w');
  fprintf(fid, '%d %d %d %d %.12f %.12f %.14f\n', dat');
  fclose(fid);
end
