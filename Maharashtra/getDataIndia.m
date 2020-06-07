function [country,C,date0] = getDataIndia()
%GETDATAINDIA Coronavirus data for India
%  as reported by One World in Data
%     https://ourworldindata.org/coronavirus-source-data
country = '(State: Maharashtra, Country: India)';
C = [
14
32
38
41
44
48
52
64
74
97
107
122
125
156
186
203
220
302
335
423
487
635
747
868
1018
1135
1364
1574
1761
1982
2334
2680
2916
3201
3321
3648
4200
4666
5218
5649
6427
6817
7628
8068
8590
9318
9915
10498
11506
12296
12974
14541
15525
16758
17974
19063
20228
22171
23401
24427
25922
27524
29100
30706
33053
35058
37136
39297
41642
44582
47190
50231


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
