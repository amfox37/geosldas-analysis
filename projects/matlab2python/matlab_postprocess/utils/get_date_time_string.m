
function [date_time_string] = get_date_time_string(date_time, file_tag) 

% reichle, 29 Jun 2005
% de lannoy, 29 Jun 2010 - updated for rst 

date_time_string = num2str(date_time.year,  '%4.4d');

if ~isempty(findstr(file_tag,'pentad')) 
  
  date_time_string = [date_time_string, 'p',                ...
		      num2str(date_time.pentad, '%2.2d')];
  
else
  
  date_time_string = [date_time_string, num2str(date_time.month, '%2.2d')];
  
end

if ~isempty(findstr(file_tag,'daily')) 
  
  date_time_string = [date_time_string, num2str(date_time.day,   '%2.2d')];
  
end

if ( ~isempty(findstr(file_tag,'xhourly'))   |  ...
     ~isempty(findstr(file_tag,'rst'))       |  ...
     ~isempty(findstr(file_tag,'innov'))     |  ...
     ~isempty(findstr(file_tag,'incr'))      |  ... %GDL added 19jul12
     ~isempty(findstr(file_tag,'OminusF'))   |  ... %GDL added 10nov10
     ~isempty(findstr(file_tag,'1hto3h'))    |  ... %GDL added 14may14
     ~isempty(findstr(file_tag,'ObsFcstAna'))|  ...
     ~isempty(findstr(file_tag,'inst'))       ) 
     
  date_time_string =                         ...
      [ date_time_string,                    ... 
	num2str(date_time.day,   '%2.2d'),   ...
	'_',                                 ...
	num2str(date_time.hour,  '%2.2d'),   ...
	num2str(date_time.min,   '%2.2d')  ];
  
end

% ============== EOF ==============================================
