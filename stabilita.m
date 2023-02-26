clear all
close all
clc

T = readtable('data_mekkaPodloz.txt'); 
T(1:19,:)=[];
T(:,1:9)=[];

array= table2array(T)
fvz =200;

[b_filter,a_filter] = butter(4,10/(fvz/2), 'low'); %koeeficity b, a pomocí butter filtru

copx = filtfilt(b_filter, a_filter, array(:,1));
copy = filtfilt(b_filter, a_filter, array(:,2));

rozsahx = max(array(:,1))-min(array(:,1));
rozsahy = max(array(:,2))-min(array(:,2));

prumer_mx=(rozsahx+rozsahx1)/2;
prumer_my=(rozsahy+rozsahy1)/2;


