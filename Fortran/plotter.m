close all;
clear all;
clf;

plmmmhll=dlmread('slug_sod_plm_10027.dat');
%plmvlhll=dlmread('slug_sod128_plm_vl_hll_10027.dat');
%plmmchll=dlmread('slug_sod128_plm_mc_hll_10027.dat');

%plmmmroe=dlmread('slug_sod128_plm_mm_roe_10027.dat');
%plmvlroe=dlmread('slug_sod128_plm_vl_roe_10027.dat');
%plmmcroe=dlmread('slug_sod128_plm_mc_roe_10027.dat');

figure(1);
hold on;
plot(plmmmhll(:,1),plmmmhll(:,2),'r');
%plot(plmvlhll(:,1),plmvlhll(:,2),'g');
%plot(plmmchll(:,1),plmmchll(:,2),'b');


%figure(2);
%hold on;
%plot(plmmmroe(:,1),plmmmroe(:,2),'r');
%plot(plmvlroe(:,1),plmvlroe(:,2),'g');
%plot(plmmcroe(:,1),plmmcroe(:,2),'b');