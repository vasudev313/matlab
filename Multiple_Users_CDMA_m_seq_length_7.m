%*************DS-CDMA Multiple User using BPSK*****************************
%*************The users considered is 5 users******************************
%*************The number of users should above ONE*************************
clear;
clc;
close all;  
Lc=7; % Spreading code length is 7
Ldata=1000000;
user=5;
BER1=[];
BER2=[];
BER_anat=[];
%*************Data Generation of Multiple User**************
data_sym=[2*round(rand(Ldata,user))-1]';
%*************m-sequence spreading code generation**************
hpn = comm.PNSequence('Polynomial',[3 2 0],'SamplesPerFrame', 7, 'InitialConditions',[1 0 0]);
pnseq= step(hpn);
pncode=(2*(pnseq)-1)';
%*********m-sequence code cyclic shifting*********
pncode_shift=zeros(user,length(pncode));
for ishift=1:user
pncode_shift(ishift,:)=(1/sqrt(Lc))*[circshift(pncode,[0,ishift])];
ishift=ishift+1;
end
%*************Spreading process***************
spread=zeros(user,Ldata*Lc);
for ispread=1:user
spread(ispread,:)=kron(data_sym(ispread,:),pncode_shift(ispread,:));
ispread=ispread+1;
end

combined_data=sum(spread);
%************Noise Generation***************
noiseq=randn(1,length(combined_data));

for iter=1:13,         
    Eb2N(iter)=(iter-1);
    Eb2N_num=10^(Eb2N(iter)/10);   
    Var_n=1/(2*Eb2N_num);      
    signois=sqrt(Var_n);        
    awgnois=signois*noiseq;     
    y_out=combined_data+awgnois;
    yreshape=reshape(y_out,Lc,Ldata); 
%*******Despreading with code***************
    z_out=pncode_shift*yreshape;
%*******Demodulating and Comparing received data with user data*************    
    dec=sign(z_out);
    todec=(dec+1)/2;
    toLdata=(data_sym+1)/2;
    
    for iber=1:user
    [number(iber),ratio(iber)]=biterr(toLdata(iber,:),todec(iber,:));
    iber=iber+1;
    end    
    BER1=[BER1;ratio];
        
%***********Q analytical********************
    BER_anat=[BER_anat;0.5*erfc(sqrt(Eb2N_num))];
end
%************To Plot************************
for ifig=1:user
    figure(ifig)
    figure_ber=semilogy(Eb2N,BER_anat,'k--',Eb2N,BER1(:,user),'r-o');
    legend('BER Analytical Curve',sprintf('User- %g',ifig));
    axis([0 12  0.99e-5 1.e0]);
    set(figure_ber,'LineWidth',2.5);
    xlabel('E_b/N (dB)'); ylabel('BPSK BER');
    title(sprintf('BER performance of User- %g',ifig))
end