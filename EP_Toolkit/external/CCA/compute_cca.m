% function [reconstr]= compute_cca(data,opt)

%written by Maarten De Vos and kindly provided to Joseph Dien 9/23/16
%included with EP Toolkit by permission of author
%10/4/17 modified to allow for specification of opt fields in function call

function [reconstr]= compute_cca(data,opt)

[a,b,c]=size(data); % dimensions of data

if nargin < 1
    opt.femg=15; % cutoff frequency of EEG between brain and muscle
    opt.fs=512; % sampling freq
    opt.range=[0 floor(a/2)]; %default component range --> never remove more than half the sources
    opt.ratio=9; % magic number :D
end;

for k=1:c  %for all the trials compute CCA
   [W(:,:,k),r] = bsscca(squeeze(data(:,:,k)),1); %effective CCA computation : W is demixing matrix
   sources(:,:,k)=squeeze(W(:,:,k))*squeeze(data(:,:,k)); %reconstructing the (CCA) sources
end



for k=1:c
    [index] = emg_psd(squeeze(sources(:,:,k)),opt) % computation of spectra and indeces of components that will be removed
    for l=1:length(index)
       sources(index(l),:,k)=zeros( 1,size(sources,2)); % zero the components that are muscle related
    end
end

for k=1:c
    M(:,:,k)=inv(squeeze(W(:,:,k))); %M is mixing matrix
    reconstr(:,:,k)=squeeze(M(:,:,k))*squeeze(sources(:,:,k)); % reconstruct the data with changed source matrix
end