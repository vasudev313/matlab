%*************DS-CDMA Single User using BPSK**************
clear;
clc;
close all;
Ldata=1000000;    
Lc=15;
BER1=[];
BER_anat=[];
%*************Data Generation of Single User**************
data_sym=[2*round(rand(Ldata,1))-1]';
%*************m-sequence spreading code generation**************
hpn = comm.PNSequence('Polynomial',[4 3 0],'SamplesPerFrame', 15, 'InitialConditions',[1 0 0 0]);
pnseq= step(hpn)';
pcode=(2*pnseq-1)/sqrt(Lc);

%*************Spreading process***************
x_in=kron(data_sym,pcode);

%************AWGN Noise Generation***************
noiseq=randn(1,Ldata*Lc); 

for iter=1:21,
    iter
    Eb2N(iter)=(iter-1);
    Eb2N_num=10^(Eb2N(iter)/10);   
    Var_n=1/(2*Eb2N_num);      
    signois=sqrt(Var_n);        
    awgnois=signois*noiseq;     

    y_out=x_in+awgnois;
    Y_out=reshape(y_out,Lc,Ldata); 
%*******Despreading with code***************
    z_out=pcode*Y_out;
    dec=sign(z_out);
%*******Demodulating received data and user data*************    
    todec=(dec+1)/2;
    toLdata=(data_sym+1)/2;
%***********Calculating BER****************
    [number,ratio]=biterr(toLdata,todec);
    BER1=[BER1;ratio];
%***********Q analytical********************
    BER_anat=[BER_anat;(1/2)*erfc(sqrt(Eb2N_num))];
end

figure(1)
figure_ber=semilogy(Eb2N,BER_anat,'k--',Eb2N,BER1,'r-o');
legend('BER Analytical Curve','User-1');
axis([0 12 0.99e-5 1.e0]);
set(figure_ber,'LineWidth',2.5);
xlabel('E_b/N (dB)'); ylabel('BPSK BER');
title('BER performance of Single-user CDMA')
    