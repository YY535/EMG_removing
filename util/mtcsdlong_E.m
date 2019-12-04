function [y2, f, t, phi, FStats, y2norm] = mtcsdlong_E(varargin);
%Modified version of mtchglong.m, which returns cross-spectrogram y2 and a product of auto-spectrograms y2norm, 
%instead of coherogram. The coherogram can be computed then as coh = y2./sqrt(y2norm).
%These cross-spectra and autospectra can be used to compute mean coherence:
% mean_coh = mean(cross-spectra) / mean(product of auto-spectra).
%
%function [yo, fo, to, phi, FStats]=mtchglong(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,FreqRange);
% Multitaper Time-Frequency Spectrum/Coherence
% for long files - splits data into blockes to save memory
% function A=mtcsg(x,nFFT,Fs,WinLength,nOverlap,NW,nTapers)
% x : input time series
% nFFT = number of points of FFT to calculate (default 1024)
% Fs = sampling frequency (default 1250)
% WinLength = length of moving window (default is nFFT)
% nOverlap = overlap between successive windows (default is WinLength/2)
% NW = time bandwidth parameter (e.g. 3 or 4), default 3
% nTapers = number of data tapers kept, default 2*NW -1
% Detrend = method of detrending before FFR (default is linear)
% output yo is yo(f, t)
%
% If x is a multicolumn matrix, each column will be treated as a time
% series and you'll get a matrix of coherence out yo(f, t, Ch1, Ch2)
% Original code by Partha Mitra - modified by Ken Harris 
% and adopted for long files and phase by Anton Sirota

% default arguments and that
[x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,nChannels,nSamples,nFFTChunks,winstep,select,nFreqBins,f,t] = mtparam(varargin);

% allocate memory now to avoid nasty surprises later
y2 = complex(zeros(nFFTChunks,nFreqBins, nChannels, nChannels)); % output array
%y=complex(zeros(nFFTChunks,nFreqBins, nChannels, nChannels)); % output array

if nargout>3
    phi=complex(zeros(nFFTChunks,nFreqBins, nChannels, nChannels));
end

%NEW output parameters: cross-spectra and normalizing factor
if nargout>5
    y2norm = complex(zeros(nFFTChunks,nFreqBins, nChannels, nChannels));
end



nFFTChunksall= nFFTChunks;
%freemem = FreeMemory;
BlockSize = 2^8;
nBlocks = ceil(nFFTChunksall/BlockSize);

% textprogressbar(['Computing cross-spectrogram: ']);

%h = waitbar(0,'Wait..');
for Block=1:nBlocks
    %   waitbar(Block/nBlocks,h);
%     textprogressbar(100*Block/nBlocks);
    
    minChunk = 1+(Block-1)*BlockSize;
    maxChunk = min(Block*BlockSize,nFFTChunksall);
    nFFTChunks = maxChunk - minChunk+1;
    iChunks = [minChunk:maxChunk];
    Periodogram = complex(zeros(nFreqBins, nTapers, nChannels, nFFTChunks)); % intermediate FFTs
    Temp1 = complex(zeros(nFreqBins, nTapers, nFFTChunks));
    Temp2 = complex(zeros(nFreqBins, nTapers, nFFTChunks));
    Temp3 = complex(zeros(nFreqBins, nTapers, nFFTChunks));
    eJ = complex(zeros(nFreqBins, nFFTChunks));
    tmpy =complex(zeros(nFreqBins,nFFTChunks, nChannels, nChannels));
    % calculate Slepian sequences.  Tapers is a matrix of size [WinLength, nTapers]
    [Tapers V]=dpss(WinLength,NW,nTapers, 'calc');
    % New super duper vectorized alogirthm
    % compute tapered periodogram with FFT 
    % This involves lots of wrangling with multidimensional arrays.
    
    TaperingArray = repmat(Tapers, [1 1 nChannels]);
    for j=1:nFFTChunks
        jcur = iChunks(j);
        if length(nOverlap)==1
        	Seg = [(jcur-1)*winstep+1: min((jcur-1)*winstep+WinLength,nSamples)];
            Segment = x(Seg,:);
        else
        	Seg = [nOverlap(jcur)-WinLength/2+1: min(nOverlap(jcur)+WinLength/2,nSamples)];
            Segment = x(Seg,:);
        end
                
        if (~isempty(Detrend))
            Segment = detrend(Segment, Detrend);
        end;
        SegmentsArray = permute(repmat(Segment, [1 1 nTapers]), [1 3 2]);
        TaperedSegments = TaperingArray .* SegmentsArray;
        
        fftOut = fft(TaperedSegments,nFFT);
        normfac = sqrt(2/nFFT); %to get back rms of original units
        Periodogram(:,:,:,j) = fftOut(select,:,:)*normfac; 
        % Periodogram: size  = nFreqBins, nTapers, nChannels, nFFTChunks
    end	
    if nargout>4 %ccompute fstats
        U0 = repmat(sum(Tapers(:,1:2:end)),[nFreqBins,1,nChannels,   nFFTChunks]);
        Mu = sq(sum(Periodogram(:,1:2:end,:,:) .* conj(U0), 2) ./  sum(abs(U0).^2, 2));
        Num = abs(Mu).^2;
        Sp = sq(sum(abs(Periodogram).^2,2));
        chunkFS = (nTapers-1) * Num ./ (Sp ./ sq(sum(abs(U0).^2, 2))- Num );
        %	sum(abs(Periodogram - U0.*repmat(Mu,[1,nTapers,1,1])), 2);
        FStats(iChunks, :, :)  = permute(reshape(chunkFS, [nFreqBins, nChannels, nFFTChunks]),[ 3 1, 2]);
    end
    % Now make cross-products of them to fill cross-spectrum matrix
    for Ch1 = 1:nChannels
        for Ch2 = Ch1:nChannels % don't compute cross-spectra twice
            Temp1 = reshape(Periodogram(:,:,Ch1,:), [nFreqBins,nTapers,nFFTChunks]);
            Temp2 = reshape(Periodogram(:,:,Ch2,:), [nFreqBins,nTapers,nFFTChunks]);
            Temp2 = conj(Temp2);
            Temp3 = Temp1 .* Temp2;
            eJ=sum(Temp3, 2);
            tmpy(:,:, Ch1, Ch2)= eJ/nTapers;
            
            % for off-diagonal elements copy into bottom half of matrix
            if (Ch1 ~= Ch2)
                tmpy(:,:, Ch2, Ch1) = conj(eJ) / nTapers;
            end            
            
        end
    end
    
    for Ch1 = 1:nChannels
        for Ch2 = 1:nChannels % don't compute cross-spectra twice
            
            if (Ch1 == Ch2)
                % for diagonal elements (i.e. power spectra) leave unchanged
                y2(iChunks,:,Ch1, Ch2) = permute(abs(tmpy(:,:,Ch1, Ch2)),[2 1 3 4]);
                
            else
                       
                % for off-diagonal elements, scale
                %y(iChunks,:,Ch1, Ch2) = permute( abs(tmpy(:,:,Ch1, Ch2)) ./ sqrt(tmpy(:,:,Ch1,Ch1) .* tmpy(:,:,Ch2,Ch2)), [2 1 3 4]);                
                phi(iChunks,:,Ch1,Ch2) = permute( angle(tmpy(:,:,Ch1, Ch2)), [2 1 3 4]);      
                
                %NEW OUTPUT PARAMETR: cross-spectra and norm. factor
                y2(iChunks,:,Ch1, Ch2) =  permute( abs(tmpy(:,:,Ch1, Ch2)) , [2 1 3 4]);
                y2norm(iChunks,:,Ch1, Ch2) = permute( tmpy(:,:,Ch1,Ch1) .* tmpy(:,:,Ch2,Ch2) , [2 1 3 4]);  
                      
            end
        end
    end %loop across Ch1
    
    
end %loop across Blocks
% textprogressbar(' DONE.');

%close(h);
% we've now done the computation.  the rest of this code is stolen from
% specgram and just deals with the output stage

if nargout == 0
    % take abs, and use image to display results
    newplot;
    for Ch1=1:nChannels, for Ch2 = Ch1:nChannels
            subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
            if Ch1==Ch2
                if length(t)==1
                    imagesc([0 1/f(2)],f,20*log10(abs(y2(:,:,Ch1,Ch2))+eps)');axis xy; colormap(jet);
                else
                    imagesc(t+diff(t(1:2))/2,f,20*log10(abs(y2(:,:,Ch1,Ch2))+eps)');axis xy; colormap(jet);colorbar
                end
                title(['Power specgram ' num2str(Ch1)]);
            else
                %imagesc the coherogram
                imagesc(t+diff(t(1:2))/2,f,(abs(y2(:,:,Ch1,Ch2)))');axis xy; colormap(jet);colorbar
                title(['Coherogram ' num2str(Ch1)]);
                
                
                %display phaseogram
                subplot(nChannels, nChannels, Ch2 + (Ch1-1)*nChannels);
                imagesc(t+diff(t(1:2))/2,f,squeeze(phi(:,:,Ch1,Ch2))');axis xy; colormap(jet);colorbar
                title(['Phasogram ' num2str(Ch1)]);
            end
        end; end;
    xlabel('Time')
    ylabel('Frequency')
end
