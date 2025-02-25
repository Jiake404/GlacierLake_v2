function [clip_hour1,clip_hour2,clip_day1,clip_day2,run_day,run_hour] = readtime(AWS,run_start,run_end)

t = readtable([AWS '_hour.csv'],'Range','A1:A2');
t = table2array(t);
AWS_start = [t.Year,t.Month,t.Day,t.Hour];
AWS_start_doy = day(datetime(AWS_start(1:3)),'dayofyear');
run_start_doy = day(datetime(run_start),"dayofyear");
run_end_doy = day(datetime(run_end),"dayofyear");
% AWS起始时间和运行起始时间之间的天数，前面需要砍去的天数
c_year = AWS_start(1):run_start(1);
if length(c_year)==1
    clip_day = run_start_doy-AWS_start_doy;
elseif length(c_year)==2
    clip_day = yeardays(AWS_start(1))-AWS_start_doy+run_start_doy;
elseif length(c_year)>2
    c_day = sum(yeardays(c_year(2:length(c_year)-1)));
    clip_day = yeardays(AWS_start(1))-AWS_start_doy+c_day+run_start_doy;
end
% AWS小时分辨率 需要砍去的小时数 确保砍去的时候从0点开始
clip_hour1 = clip_day*24-AWS_start(4)+1;
% 运行天数
run_day = length(datetime(run_start):datetime(run_end));
% 运行小时数
run_hour = run_day*24;
clip_hour2 = clip_hour1+run_hour-1;

%RACMO(从1-1起)裁剪
% 计算起始时间之间有多少天
num_day = datenum(run_end) - datenum(run_start);
clip_day1 = run_start_doy;
clip_day2 = run_start_doy+num_day;

