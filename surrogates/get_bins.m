function bins=get_bins(filename,dim,lag);

% function bins=get_bins(filename,dim,lag);
%
% Opensthe file [filename,'.',dim,'.',lag,'.bin] looks for the line
% starting with '# scale factor used:' and then returns the rest of the line 
% as bins. is bins is the scale factor used in the binning produced by
% the file  [filename,'.',dim,'.',lag,'.bin].
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk


bins=[];
if lag>1,
  filename=[filename,'.',int2str(dim),'.',int2str(lag),'.bin'];
else
  filename=[filename,'.',int2str(dim),'.bin'];
end;
[fid,message]=fopen(filename);
if fid~=-1
	i=0;
	while 1
		line=fgetl(fid);
		if ~isstr(line), break, end
		if length(line)>19
		if strcmp(line(1:20),'# scale factor used:')
			disp(line);
			bins=eval(line(21:length(line)));
			fclose(fid);
			return
			disp('This shouldnt happen');
		end;
		end;
	end;
	fclose(fid);
else
	disp(message);
end;





