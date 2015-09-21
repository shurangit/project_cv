function [video_path,img_files] = load_video(base_path)
%CHOOSE_VIDEO

	if ispc(), base_path = strrep(base_path, '\', '/'); end
	if base_path(end) ~= '/', base_path(end+1) = '/'; end
	
	%list all sub-folders
	contents = dir(base_path);
	names = {};
	for k = 1:numel(contents),
		name = contents(k).name;
		if isdir([base_path name]) && ~strcmp(name, '.') && ~strcmp(name, '..'),
			names{end+1} = name;  %#ok
		end
	end
	
	%no sub-folders found
	if isempty(names), video_path = []; return; end
	
	%choice GUI
	choice = listdlg('ListString',names, 'Name','Choose video', 'SelectionMode','single');
	
	if isempty(choice),  %user cancelled
		video_path = [];
	else
		video_path = [base_path names{choice} '/'];
	end
	


text_files = dir([video_path '*_frames.txt']);
    if ~isempty(text_files),
		f = fopen([video_path text_files(1).name]);
		frames = textscan(f, '%f,%f');
		fclose(f);
		
		%see if they are in the 'imgs' subfolder or not
		if exist([video_path num2str(frames{1}, 'imgs/img%05i.png')], 'file'),
			video_path = [video_path 'imgs/'];
		elseif ~exist([video_path num2str(frames{1}, 'img%05i.png')], 'file'),
			error('No image files to load.')
		end
		
		%list the files
		img_files = num2str((frames{1} : frames{2})', 'img%05i.png');
		img_files = cellstr(img_files);
	else
		%no text file, just list all images
		img_files = dir([video_path '*.png']);
		if isempty(img_files),
			img_files = dir([video_path '*.jpg']);
			assert(~isempty(img_files), 'No image files to load.')
		end
		img_files = sort({img_files.name});
    end
end