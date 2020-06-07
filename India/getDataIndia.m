function [country,C,date0] = getDataIndia()
%GETDATAINDIA Coronavirus data for India
%  as reported by One World in Data
%     https://ourworldindata.org/coronavirus-source-data
country = '(Country: India)';
C = [
102
112
126
146
171
198
256
334
403
497
571
657
730
883
1019
1139
1326
1635
2059
2545
3105
3684
4293
4777
5350
5915
6728
7599
8453
9211
10454
11485
12371
13432
14354
15725
17305
18544
20081
21373
23040
24448
26283
27890
29458
31360
33065
34866
37262
39826
42778
46434
49405
53007
56351
59690
62865
67176
70767
74330
78056
82047
85855
90649
95698
100326
106480
112200
118223
124759
131422
138535


%           1 % 21-Mar-2020
%           6 % 04-Mar-2020
%          28 % 05-Mar-2020
%          29 % 06-Mar-2020
%          31 % 07-Mar-2020
%          34 % 08-Mar-2020
%         NaN % 09-Mar-2020
%          44 % 10-Mar-2020
%          50 % 11-Mar-2020
%          73 % 12-Mar-2020
%          75 % 13-Mar-2020
%          83 % 14-Mar-2020
%          90 % 15-Mar-2020
%          93 % 16-Mar-2020
%         125 % 17-Mar-2020
%         137 % 18-Mar-2020
%         165 % 19-Mar-2020
%         191 % 20-Mar-2020
%         231 % 21-Mar-2020
%         320 % 22-Mar-2020
%         439 % 23-Mar-2020
%         492 % 24-Mar-2020
%         562 % 25-Mar-2020
%         649 % 26-Mar-2020
%         724 % 27-Mar-2020
%         873 % 28-Mar-2020
%         979 % 29-Mar-2020
%        1071 % 30-Mar-2020
%        1251 % 31-Mar-2020
%        1397 % 01-Apr-2020
%        1965 % 02-Apr-2020
%        2301 % 03-Apr-2020
%        2902 % 04-Apr-2020
%        3374 % 05-Apr-2020
%        4067 % 06-Apr-2020
%        4421 % 07-Apr-2020
%        5194 % 08-Apr-2020
%        5734 % 09-Apr-2020
%        6412 % 10-Apr-2020
%        7447 % 11-Apr-2020
%        8356 % 12-Apr-2020
%        9152 % 13-Apr-2020
%       10363 % 14-Apr-2020
%       11438 % 15-Apr-2020
%       12380 % 16-Apr-2020
%       13387 % 17-Apr-2020
%       14378 % 18-Apr-2020
%       15712 % 19-Apr-2020
%       17265 % 20-Apr-2020
%       18600 % 21-Apr-2020
%       19984 % 22-Apr-2020
%       21393 % 23-Apr-2020
%       23077 % 24-Apr-2020
%       24506 % 25-Apr-2020
%       26496 % 26-Apr-2020
%       27892 % 27-Apr-2020
%       29435 % 28-Apr-2020
%       30892 % 29-Apr-2020
     % 3108 % 27-Apr-2020
%<-------------- add new data here
]';
date0=datenum('14-Mar-2020');
end
