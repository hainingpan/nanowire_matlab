% [dosmap,rev,dosmat]=spec_seinhom_sp(1,1,.2,5,.2,inf,300,'const',0,0,0,1);
% randlist=load('D:\CMTC\bothlead\Rp25\disorder\1\muVarlist.dat');
% spec_sedis_sp(1,1,.2,5,.2,inf,300,0.2,randlist,1);
% randlist=load('D:\CMTC\bothlead\Rp25\good\muVar\L3\randlist.dat');
% spec_sedis_sp(1,1,.2,5,.2,inf,300,0.4,randlist,0);
% spec_sedis_sp(1,1,.2,5,.2,inf,300,0.4,randlist,1);
% randlist=load('D:\CMTC\bothlead\Rp25\disorder\61\muVarlist.dat');
% spec_sedis_sp(1,1,.2,5,.2,inf,300,0.6,randlist,1);
% randlist=load('D:\CMTC\bothlead\Rp25\disorder\23\muVarlist.dat');
% spec_sedis_sp(1,1,.2,5,.2,inf,300,0.8,randlist,1);
% randlist=load('D:\CMTC\bothlead\Rp25\ugly\muVar\L3\randlist.dat');
% spec_sedis_sp(1,1,.2,5,.2,inf,300,1,randlist,0);
% spec_sedis_sp(1,1,.2,5,.2,inf,300,1,randlist,1);
% randlist=load('D:\CMTC\bothlead\Rp25\disorder\4\muVarlist.dat');
% spec_sedis_sp(1,1,.2,5,.2,inf,300,2,randlist,0);
% spec_sedis_sp(1,1,.2,5,.2,inf,300,2,randlist,1);

% 
% randlist=load('D:\CMTC\bothlead\Rp25\good\DeltaVar\L3\randlist.dat');
% spec_segap_sp(1,1,.2,5,.2,inf,300,.06,randlist,1);
% randlist=load('D:\CMTC\bothlead\Rp28\gVar\3\randlist.dat');
% spec_seg_sp(1,1,.2,5,.2,inf,300,0.3,randlist,1);
% spec_seqd_sp(1,1,.2,5,.2,inf,0.6,40,300,1);
% spec_seinhom_sp(1,1,.2,5,.2,inf,300,'exp',-1.2,0,40,1);
% randlist=load('D:\CMTC\bothlead\Rp28\QD+muVar\1\muVarlist0.1.dat');
% [dosmap,rev,dosmat]=spec_seqddis_sp(1,1,.2,5,.2,inf,.6,40,300,.1,randlist,1);
% randlist=load('D:\CMTC\bothlead\Rp28\QD+muVar\1\muVarlist1.dat');
% spec_seqddis_sp(1,1,.2,5,.2,inf,.6,40,300,1,randlist,1);
% randlist=load('D:\CMTC\bothlead\Rp28\QD+muVar\23\muVarlist2.dat');
% spec_seqddis_sp(1,1,.2,5,.2,inf,.6,40,300,2,randlist,1);

% randlist=load('d:\cmtc\bothlead\Rp17\ugly\1\randlist.dat');
% spec_sedis_sp(1,1,.2,5,.2,1.2,100,1,randlist);
% randlist=load('d:\cmtc\bothlead\Rp17\ugly\3\randlist.dat');
% spec_sedis_sp(1,1,.2,5,.2,1.2,100,1,randlist);
% randlist=load('d:\cmtc\bothlead\Rp16\ugly\1\vimp.txt');
% spec_sedis_sp(1,1,.2,5,.2,1.2,300,1,randlist);
% randlist=load('d:\cmtc\bothlead\Rp19\gvar1\104\randlist.dat');
% spec_seg_sp(1,1,.2,5,.2,1.2,100,0.8,randlist);
% randlist=load('d:\cmtc\bothlead\Rp19\gvar1\115\randlist.dat');
% spec_seg_sp(1,1,.2,5,.2,1.2,100,0.8,randlist);
% randlist=load('d:\cmtc\bothlead\Rp19\gvar3\13\vimp.txt');
% spec_seg_sp(1,1,.2,5,.2,1.2,300,0.6,randlist);

for i=2:10
    [rev,re,vimp]=spec_dis(1,i,300,0,0);
end
