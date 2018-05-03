close all
clear all
clc

%% reading the original song and downsampling to 1000
[song,F]=audioread('test_audio.wav')  ;    %% reading the audio file
info_1 = audioinfo('test_audio.wav')  ;    %% get information of the audio file
 sound (song,F) ;    %% play the audio file
figure
plot(song) ;   %%plot the original signal
title('Original Signal input');
resam=resample(song,F,((F/1000)*F)); %% downsampling to 1000
figure
plot(resam); %% plot the downsampled signal
title('Down Sampled Signal  input');
wavwrite(resam,1000,'song_resamplees.wav');   %% write the downsampled signal as a wav file
[re,f]=audioread('song_resamplees.wav');      %%  load the resampled audio file
info_2 = audioinfo('song_resamplees.wav');    %% get information of the downsampled audio file
% sound (re,f)     %% play the  downsampled audio file

%% Quantization Process
Maximum_value_for_the_signal=max(max(resam)'); %% maximum value for the signal
Minimum_value_for_the_signal=min(min(resam)'); %% minimum value for the signal
resam=resam-Minimum_value_for_the_signal; %% signal is shifted above the zero to make it easier for mapping with levels                                                                        %number of bits of quantizer
delta=(Maximum_value_for_the_signal-Minimum_value_for_the_signal)/((2^13)-1);  %% the step size between each two levels
Analog_levels = (0:8191)'; %% look up table for saving the binary corresponding to each level value
Digital_levels = de2bi(Analog_levels); %% corresponding binary to each level
New= [] ;
New_matrix_for_saving_binary_data=[ ] ; %% empty matrix for saving the binary data
Matrix_for_first_520_data=[];
IIi=1;
index=1;

for i=1:length(resam) %% loop for every element on the resamples signal
    for Index=1:(2^13) %% loop for every level from the lookup table
        if(resam(i,1)>=(Index*delta))&&(resam(i,1)<=((Index+1)*delta)) %% find the two levels that include the value between them
            if(resam(i,1)<((Index*delta)+(delta/2))) %% if the resapled value is less than the threshold between the two vales
                New_matrix_for_saving_binary_data(i,:)=Digital_levels(Index,:); %% then approximate it to the nearst lower level
                break;
            else %% if the resapled value is greater than the threshold between the two vales
                New_matrix_for_saving_binary_data(i,:)=Digital_levels((Index+1),:); %% then approximate it to the nearst higher level
                break;
            end
        else %% loop again and take another two levels to compare
        end
    end
end
[rows,cols]=size(New_matrix_for_saving_binary_data); %% finding the size of the generated matrix

% concatenating all the data in one vector
Vector_containing_all_data=[]; %% empty matrix for concatenating the data
for Io=1:rows
    Vector_containing_all_data=[Vector_containing_all_data New_matrix_for_saving_binary_data(Io,:)];
end
% to make it possible for the generated vector to be divided to 260's bits
% we have to add some pooling data
Polling_zeros=zeros(1,1157);
Vector_containing_all_data_after_adding_polling=[Vector_containing_all_data Polling_zeros ]; %%concatenating the pooling zeros to the data
[row_div,col_div]=size(Vector_containing_all_data_after_adding_polling);
%%

Number_of_260_bits=col_div/260; %% after dividing the data to multiple vectors of 260 bits
Used_260=[]; %% empty matrix to save the used 260 bit vector
z=0; %% intial value for calculating the 260 bits for each vector
B=[]; %% empty matrix for saving all frames
concatentation_demodulated_signal=[]; %% empty matrix for concatenating the received 156 bits in one vector
All_Mod=[ ] ; %% matrix for concatenating the modulated data
All_Noisy_Mod=[ ] ; %% matrix for concatenating the modulated data with adding the noise effect
Channel_Coding=[ ] ;
for Ii=1:13 %% loop for the number of second
    for Ii1=1:25 %% loop for number of milli seconds in one second = 50 we have divided to 2 for loops to work for eash 25 milii second separately
        for Ii2=1:2
            %% channel encoder
            QQ=1; %% intial value for the index to take each 260 bit separately
            Used_260=Vector_containing_all_data_after_adding_polling((QQ+z):(QQ+259+z)) ;
            z=z +260 ;
            c1a =Used_260(1:50); %% class 1 a containing the first 50 bits of data
            c1b =Used_260(51:182); %% class 1 b containing the next 130 bits of data
            c2 =Used_260(183:260); %% class 2 containing the next 77 bits of data
            %PARITY ENCODING:THREE CHECK BITS ARE ADDED
            g = [1 0 1 1];
            gen = comm.CRCGenerator('Polynomial',g,'ChecksumsPerFrame',1);
            codeword = step(gen,c1a'); %% three check bits are generated
            c1 = [codeword' c1b 0 0 0 0]; %% concatenating class 1a , class 2a and four tail bits are added finally
            constraint_length=5; %% number of registers
            polynomial=[25 35];
            trellis=poly2trellis(constraint_length,polynomial); %% convert to trellis diagram using the generator polynomial
            enc_c1=convenc(c1,trellis);
            ch_coding_data = [enc_c1 c2]; %% concatenating the encoded bits and class 
            Channel_Coding=[ Channel_Coding ch_coding_data];
        
            %% First Level Interleaving 
            
            Frames=[ ] ; %% empty matrix to save the frames
            u=1; %% intial index for saving each 8 bits on a frame
            for l=1:57
                frame_formation=ch_coding_data(u :u +7 ) ; %% concatenating each 8 bits
                Frames=[ Frames ; frame_formation ] ; %% forming a frame matrix of size 57*8
                u= u + 8 ;
            end
            B=[B Frames ]; %% concatenating all frames in one matrix
        end
        % Burst Format and Second Level Interleaving 
        c=1;
        e=0 ;
        Traning_bits=[1 1 1 1 1 0 0 1 1 0 1 0 1 0 0 0 0 0 1 1 0 0 1 0 1 0 ]; %% given training bits
        Tails_bits=[0 0 0]; %% three zeors tail bits
        flag=[0] ; %% flag =0
        Burst_Format=[] ;
        Gard_bits=zeros(1,8); %% eight quard zero bits
        
        for i=1:8
            burstframe_1=B(:,i); %% taking the first column of data
            burstframe_2=B (:,i+8); %% taking the coresponding data from the next frame
            Burst=[ Tails_bits burstframe_1'  flag Traning_bits  flag  burstframe_2'  Tails_bits Gard_bits] ; %% concatenating the burst data
            
            % Modulation
            Tb=1/(270*10^3); %% Tb=1/Rb
            Fc=(270*10^3); %% Fc=270KHZ
            hMod = comm.GMSKModulator('BitInput', true, 'InitialPhaseOffset', pi/4,'SamplesPerSymbol',10);
            hAWGN = comm.AWGNChannel('NoiseMethod', ...
                'Signal to noise ratio (SNR)','SNR',25);
            Burst1=[Burst'; 0]; %% concatenating zero to the burst to match the size
            modSignal = step(hMod, Burst1);
            Length_of_t=length(modSignal);
            t = linspace(0,(Length_of_t-1)*Tb/10,Length_of_t);
            modSignal_carr = real(modSignal.*exp(1i*2*pi*Fc*t'));
            All_Mod=[ All_Mod modSignal_carr ] ;  %%% Matrix contain all modulated signal without noise 
            noisySignal = step(hAWGN, modSignal_carr);  
            All_Noisy_Mod=[ All_Noisy_Mod noisySignal ]; %%% Matrix contain all modulated signal with noise 
            
        end
        B=[ ] ; %% intiallize the matrix B for receiving new data
        Frames=[ ] ; %% intiallize the matrix Frames for receiving new data
    end
    
end
    figure 
   plot( Channel_Coding)
       title('Channel Coded Data ');
%%%%%%%%%%%%%% Reciver %%%%%%%%%%%%%%%
[XX,YY]=size(All_Mod); %% size of the received modulated matrix
for R=1:YY %% loop for all the columns 
    Mod_Signal=All_Mod(:,R);
    noisySignal=All_Noisy_Mod(:, R);
     DemodSignal_carr =  Mod_Signal.*exp(-1i*2*pi*Fc*t'); %% canceling the effect of the exponential
%        DemodSignal_carr =  noisySignal.*exp(-1i*2*pi*Fc*t'); %%% Demodulated signal with noise
    hDemod = comm.GMSKDemodulator('BitOutput', true, ...
        'InitialPhaseOffset', pi/4,'TracebackDepth',1,'SamplesPerSymbol',10);
    receivedData = step(hDemod, DemodSignal_carr);
    receivedData=receivedData((2:end),1);
    [mm,nn]=size(receivedData); %% receiving the 156 bits
    concatentation_demodulated_signal=[concatentation_demodulated_signal ; receivedData']; %% concatenating all data in one vector
end

Recived_Blocks=[ ] ; %% empty matrix for savine the data from the channel decoding
Q=0; %% intialize the index for
for b=1:13 %% loop for the number of second
    % intialize the matrices for saving data for each second
    First_Frame=[];
    Second_Frame=[ ] ;
    F1=[ ] ;
    F2=[];
    for g=1:25 %% loop for number of milli seconds in one second = 50 we have divided to 2 for loops to work for eash 25 milii second separately
        Second_Frame=[ ] ;
        First_Frame= [ ] ;
        for f=1:8 %% loop for each column of the frame
            Burst_1= concatentation_demodulated_signal(Q+f ,:);
            frame1=Burst_1(4:60)';    %%take the Frame 1 (57 bit ) from burst format 
            frame2=Burst_1(89:145)';   %%take the Frame 2 (57 bit ) from burst format 
            First_Frame=[First_Frame  frame1 ];
            Second_Frame=[ Second_Frame frame2 ];
        end
        F1=[] ;
        F2=[];
        Q=Q+8 ;
        for a=1:57 %% vector formation
            frame_formation=First_Frame (a,:); %%concatenate all coulms of Frame in one vector 
            F1=[ F1 frame_formation];
            Frame_Formation2= Second_Frame (a,:);
            F2=[F2 Frame_Formation2];
        end
        for aa=1:2 %% loop for using the first or second frame
            Frame_used=[];
            if (aa==1)
                Frame_used=F1;
            elseif(aa==2)
                Frame_used=F2;
            else
            end
            L = length(Frame_used); %% size of frame used each time
            rx_c1 = Frame_used(1:378); %% class 1 received
            rx_c2 = Frame_used(379:L); %% class 2 received
            traceback=5*constraint_length;
            dec_data=vitdec(rx_c1,trellis,traceback,'trunc','hard'); %% decoded data after the viterebi decoder
            main_dec_data=dec_data(1:185);
            rx_c1_a=main_dec_data(1:53);
            rx_c1_after_canceling_crc_bits=rx_c1_a(1:50); %% removing the three bits generated from the channel coder
            DATA_Ib = main_dec_data(54:185);
            rx_block = [rx_c1_after_canceling_crc_bits DATA_Ib rx_c2]; %% concatenating data after the channel decoding
            Recived_Blocks=[ Recived_Blocks rx_block ] ; %% concatenating all data from the channel decoder
            
        end
        
    end
end
    figure 
            plot( Recived_Blocks)
title('Channel Decoded Data  ');
[mmn,mmnn]=size(Recived_Blocks); %% size of the received data
Recived_vector_with_out_pooling=Recived_Blocks(:,(1:(mmnn-1157))); %% removing the pooling bits added above
[name_row,name_col]=size(Recived_vector_with_out_pooling); %% size of the received data after removing the pooling bits
Matrix_containing_all_data_renormalized=[]; %% empty matrix for saving the data each 13 bits in one row to make ot easy for the dac part
for Io=1:13:name_col
    Matrix_containing_all_data_renormalized=[Matrix_containing_all_data_renormalized ; Recived_vector_with_out_pooling(:,Io:(Io+12))];
end

[Nb,Mb]=size(Matrix_containing_all_data_renormalized); %% size of matrix after arranging every 13 bits in one row
Analog_levels_again=[]; %% empty matrix to save the corresponding decimal values to each binary row
for KkK=1:Nb
    Analog_levels_again(KkK,:)=bi2de(Matrix_containing_all_data_renormalized(KkK,:))+1;
end

% finding the analog value corresponding to each level to reassign the
% analog values to each level
corresponding_analog_value_to_each_level=[]; %% empty matrix to save the analog values corresponding to each level
for Jkl=1:8192
    corresponding_analog_value_to_each_level=[corresponding_analog_value_to_each_level;(Jkl*delta)];
end

[size_row,size_col]=size(Analog_levels_again);
Signal_reconstructed=[]; %% empty matrix to assign the analod values to each element of the received data
for Yun=1:Nb
    Signal_reconstructed(Yun,:)=corresponding_analog_value_to_each_level(Analog_levels_again(Yun,:),:);
end
Signal_reconstructed=Signal_reconstructed+Minimum_value_for_the_signal; %% bring the signal back to it's position by adding the minimum value
Recieved_signal=Signal_reconstructed;

figure
plot(Recieved_signal); %% plot the downsampled signal
title('Down sampled signal output');
wavwrite(Recieved_signal,1000,'final_down_sampled.wav'); %% listening to the downsampled received signal
[ree,ff]=audioread('final_down_sampled.wav');
info_22 = audioinfo('final_down_sampled.wav');
% sound (ree,ff)

% upsampling process equivalent to the speech decoding part
Recieved_signal_1=resample(Recieved_signal,44100,(1/44.1)*44100);
figure
plot(Recieved_signal_1); %% plot the Reconstructed  signal
title('Reconstructed Signal ');

% listening to the received signal after the upsampling process
wavwrite(Recieved_signal_1,44100,'Reconstructed_signal .wav');
[reee,fff]=audioread('Reconstructed_signal .wav');
info_22 = audioinfo('Reconstructed_signal .wav');
sound (reee,fff)