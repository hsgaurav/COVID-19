close all
nIndia =26;
res = fitVirusCV19(@getDataIndia,nIndia,'prn','on','jpg','on','jpres','-r300');
res = fitVirusCV19R(@getDataIndia,'prn','on','jpg','on','jpres','-r300');

