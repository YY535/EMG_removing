function [yo, fo] = mtcsdfast(varargin);
% function [yo, fo]=mtcsdfast(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers, FreqRange, TriggTimes, ChannelPairs);
%
%A modification of the original Multitaper Cross-Spectral Density function mtcsd.m which:
% 1) exploits 4D matrix-vectorized algorithm instead of a loop over FFTChunks. It is faster than mtcsd.m, but requires more memory.
% 2) accepts as an optional input parameter a vector of triggering timestamps. If the triggering timestamps are provided, 
%     all windows will be centered at these times.
% 3) accepts as an optional input parameter ChannelPairs - a Nx2 matrix with pairs of channels for which cross-spectra must be computed.
%Evgeny Resnik, v. 03.04.2012
%
% function [yo, fo]=mtcsd(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers, FreqRange);
% Multitaper Cross-Spectral Density
% x : input time series
% nFFT = number of points of FFT to calculate (default 1024)
% Fs = sampling frequency (default 2)
% WinLength = length of moving window (default is nFFT)
% nOverlap = overlap between successive windows (default is WinLength/2)
% if nOverlap is a vector - then it's gibing times of the window centers to
% compute the spectral estimate over. units - same sampling rate as signal
% NW = time bandwidth parameter (e.g. 3 or 4), default 3
% nTapers = number of data tapers kept, default 2*NW -1
%
% output yo is yo(f)
%
% If x is a multicolumn matrix, each column will be treated as a time
% series and you'll get a matrix of cross-spectra out yo(f, Ch1, Ch2)
% NB they are cross-spectra not coherences. If you want coherences use
% mtcohere

% Original code by Partha Mitra - modified by Ken Harris
% Also containing elements from specgram.m


% default arguments and that
% [x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,nChannels,nSamples,nFFTChunks,winstep,select,nFreqBins,f,t, FreqRange] = mtparam(varargin);
[x, nFFT, Fs, WinLength, nOverlap, NW, Detrend, nTapers, nChannels, nSamples, nFFTChunks, winstep, select, nFreqBins, f, t, FreqRange, TriggTimes, ChannelPairs] ...
    = mtparam_trigg(varargin);


clear varargin; % since that was taking up most of the memory!

% check for column vector input
if nSamples == 1
	x = x';
	nSamples = size(x,1);
	nChannels = 1;
end;

if any(unique(ChannelPairs) > nChannels)  
    error('Signal has fewer channels than specified in ChannelPairs!')
end

% calculate Slepian sequences.  Tapers is a matrix of size [WinLength, nTapers]
if nFFTChunks==1
%     [Tapers V]=dpss(nSamples,NW,nTapers, 'calc');  %original, which doesn't match a size of SegmentsArray later, when nFFTChunks=1 and WinLength<nSamples
    [Tapers V]=dpss(WinLength,NW,nTapers, 'calc'); %temporal solution by Evgeny. Must be tested. ??? QUESTION What is different here?
else
    [Tapers V]=dpss(WinLength,NW,nTapers, 'calc');
end


% New super duper 4D vectorized alogirthm
% compute tapered periodogram with FFT 
% This involves lots of wrangling with multidimensional arrays.

%Index matrix for sliding windows [Time x FFTChunks]
if isempty(TriggTimes)
    %No triggering times provided. Use conventional successive sliding windows with an overlap
    FFTChunkWin = repmat([1:WinLength],nFFTChunks,1) + repmat([0:winstep:(nFFTChunks-1)*winstep]',1,WinLength);
    FFTChunkWin = FFTChunkWin';
    
else %NEW
     %Triggering times provided. Use sliding windows centered at the triggering times.
     %Ignore: nOverlap, winstep
     %Recompute: nFFTChunks
     disp('Triggering times are provided. Sliding windows will be centered at the triggering times!')
     
     %length of a sliding window halves around the triggering time  (assymetric in case of the even WinLength)
     nBefore = ceil((WinLength-1)/2);
     nAfter   = fix((WinLength-1)/2);
     %discard triggering times with a window limits outside the signal (all windows must be of the same size)
     BadTriggTimes =  TriggTimes-nBefore<1  |  TriggTimes+nAfter>nSamples;
     TriggTimes(BadTriggTimes)=[];
     clear BadTriggTimes     
          
     %recalculate the number of chunks, ignoring the output of mtparam_trigg
     nFFTChunks = length(TriggTimes);   
     
     %Inidices of samples within each window centered at the triggering times [TriggTimes-nBefore : TriggTimes+nAfter]
     FFTChunkWin = repmat(TriggTimes, 1, WinLength)  + repmat([-nBefore : nAfter]   , length(TriggTimes), 1);   
     FFTChunkWin = FFTChunkWin';     
     clear nBefore nAfter
end


% allocate memory now to avoid nasty surprises later
Segment = zeros(nChannels, WinLength, nFFTChunks);
Periodogram = complex(zeros(nFreqBins, nTapers, nChannels, nFFTChunks)); %  size(Periodogram)
Temp1 = complex(zeros(nFreqBins, nTapers, nFFTChunks));
Temp2 = complex(zeros(nFreqBins, nTapers, nFFTChunks));
Temp3 = complex(zeros(nFreqBins, nTapers, nFFTChunks));
eJ = complex(zeros(nFreqBins,nFFTChunks));
y = complex(zeros(nFreqBins, nChannels, nChannels)); % output array


%Determine which channels must be processed based on ChannelPairs (can speed up calculation)
if isempty(ChannelPairs)
    %ChannelPairs is not provided (empty): compute over all possible channel combinations without repetitions
    UsedChannels = 1:nChannels;
    ChannelPairs = unique(sort(combn(UsedChannels,2),2),'rows'); 
else
    %ChannelPairs is provided:
    UsedChannels = unique(ChannelPairs);
end
nUsedChannels = length(UsedChannels);
nChannelPairs = size(ChannelPairs,1);


%Create a matrix with data segments corresponding to the windows
%(loop over only the specified in ChannelPairs channels, other channels stay zero)
for c=1:nUsedChannels
    %extract LFP data from the channel
    tmp = x(:, UsedChannels(c) );
    SegmentChan = tmp(FFTChunkWin);     
    %detrend extracted LFP data
    if (~isempty(Detrend))
        SegmentChan = detrend( SegmentChan  , Detrend);
    end;      
    Segment( UsedChannels(c) , : , :) = SegmentChan;    
end %loop across channels
clear tmp c SegmentChan

%original version: over all channels
% for c=1:nChannels  
%     %extract LFP data from the channel
%     tmp = x(:,c);
%     SegmentChan = tmp(FFTChunkWin);     
%     %detrend extracted LFP data
%     if (~isempty(Detrend))
%         SegmentChan = detrend( SegmentChan  , Detrend);
%     end;      
%     Segment(c, : , :) = SegmentChan;    
% end %loop across channels
% clear tmp cS egmentChan


%Tapering matrix
TaperingArray = repmat(Tapers, [1 1 nChannels nFFTChunks]); %  size(TaperingArray)
SegmentsArray = repmat(Segment, [1 1 1 nTapers]) ; %  size(SegmentsArray)
SegmentsArray = permute(SegmentsArray, [2 4 1 3]);
TaperedSegments2 = TaperingArray .* SegmentsArray; %  size(TaperedSegments2)
fftOut2 = fft(TaperedSegments2,nFFT); %  size(fftOut2)
 %norm factor to get the original units (rms)
normfac = sqrt(2/nFFT);  %disregard that 0 and Fnyq should not have 2 in the normfac
% but we are not computing them anyway
Periodogram(:,:,:,:)  = fftOut2(select,:,:,:).*normfac;


% Now make cross-products of them to fill cross-spectrum matrix
%(loop only over the specified ChannelPairs, other pairs values stay zero)
for p=1:nChannelPairs
    Ch1 = ChannelPairs(p,1);
    Ch2 = ChannelPairs(p,2);
    Temp1 = squeeze(  Periodogram(:, :, Ch1, :)  );
    Temp2 = squeeze( Periodogram(:, :, Ch2, :)  );
    Temp2 = conj(Temp2);
    Temp3 = Temp1 .* Temp2;
    eJ=squeeze(sum(Temp3, 2));
    y(:,Ch1, Ch2) = mean(eJ'/nTapers,1)'  ;     %     y(:,Ch1, Ch2)= mean(eJ'/nTapers)'   ; original   and wrong!  
    y(:,Ch2, Ch1) = conj(y(:,Ch1,Ch2)) ;         %fill other half of matrix with complex conjugate
end %loop across channel pairs
clear p Ch1 Ch2 eJ Temp*


% Calculate spectra for individual channels in ChannelPairs 
%(needed only in mtchdfast.m for normalization and converting to coherence!)
for c=1: nUsedChannels
    Temp1 = squeeze(  Periodogram(:, :,  UsedChannels(c), :)  );
    Temp2 = squeeze( Periodogram(:, :,  UsedChannels(c), :)  );
    Temp2 = conj(Temp2);
    Temp3 = Temp1 .* Temp2;
    eJ=squeeze(sum(Temp3, 2));
    y(:, UsedChannels(c),  UsedChannels(c)) = mean(eJ'/nTapers,1)'  ;
end %loop across channels used in ChannelPairs
clear c  eJ Temp*


%original version: over all channels
% for Ch1 = 1:nChannels
%     for Ch2 = Ch1:nChannels % don't compute cross-spectra twice
%         Temp1 = squeeze(Periodogram(:,:,Ch1,:));
%         Temp2 = squeeze(Periodogram(:,:,Ch2,:));
%         Temp2 = conj(Temp2);
%         Temp3 = Temp1 .* Temp2;
%         eJ=squeeze(sum(Temp3, 2));
%         y(:,Ch1, Ch2)= mean(eJ'/nTapers,1)';
%     end
% end

%original version: now fill other half of matrix with complex conjugate
% for Ch1 = 1:nChannels
% 	for Ch2 = (Ch1+1):nChannels
% 		y(:, Ch2, Ch1) = conj(y(:,Ch1,Ch2))    ;
% 	end
% end

		
% we've now done the computation.  the rest of this code is stolen from
% specgram and just deals with the output stage
if nargout == 0
	% take abs, and plot results
    %newplot;
    disp('plotting....may take a while ...')
    figure
    for Ch1=1:nChannels, for Ch2 = 1:nChannels
    	subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
        plot(f,10*log10(abs(y(:,Ch1,Ch2))+eps));
		grid on;
		if(Ch1==Ch2) 
			ylabel('psd (dB)'); 
		else 
			ylabel('csd (dB)'); 
		end;
		xlabel('Frequency');
        title(['Ch' num2str(Ch1) ' vs. Ch'  num2str(Ch2)  ])
%         axis tight
	end; end;
elseif nargout == 1
    yo = y;
elseif nargout == 2
    yo = y;
    fo = f;
end


%FOR DEBUGGING: mimic overlaping windows mode with triggering times
% nSamples = size(lfp,1);
% winstep = WinLength - nOverlap;
% nFFTChunks = max(1,round(((nSamples-WinLength)/winstep)))
% FFTChunkWin = repmat([1:WinLength],nFFTChunks,1) + repmat([0:winstep:(nFFTChunks-1)*winstep]',1,WinLength);
% nBefore = ceil((WinLength-1)/2);
% nAfter   = fix((WinLength-1)/2);
% TriggTimes = FFTChunkWin(:,nBefore+1);
% %FFTChunkWin2 = repmat(TriggTimes, 1, WinLength)  + repmat([-nBefore : nAfter]   , length(TriggTimes), 1);
% ChannelPairs = unique(sort(combn(1:size(lfp,2),2),2),'rows'); %[2 3; 3 4]
% [y1, f1] = mtcsd(lfp, nfft, lfpSamplingRate, WinLength, nOverlap, [], [], [], FreqRange); 
% [y2, f2] = mtcsdfast(lfp, nfft, lfpSamplingRate, WinLength, nOverlap, [], [], [], FreqRange); 
% [y3, f3] = mtcsdfast(lfp, nfft, lfpSamplingRate, WinLength, nOverlap, [], [], [], FreqRange, TriggTimes, ChannelPairs


%-------------------------------------------------------------------------------%
%--------------------------- Supplementary functions -------------------------------%
%-------------------------------------------------------------------------------%

% helper function to do argument defaults etc for mt functions
function [x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,nChannels,nSamples,nFFTChunks,winstep,select,nFreqBins,f,t,FreqRange, TriggTimes, ChannelPairs] ...
    = mtparam_trigg(P)

nargs = length(P);

x = P{1};
if (nargs<2 | isempty(P{2})) nFFT = 1024; else nFFT = P{2}; end;
if (nargs<3 | isempty(P{3})) Fs = 1250; else Fs = P{3}; end;
if (nargs<4 | isempty(P{4})) WinLength = nFFT; else WinLength = P{4}; end;
if (nargs<5 | isempty(P{5})) nOverlap = WinLength/2; else nOverlap = P{5}; end;
if (nargs<6 | isempty(P{6})) NW = 3; else NW = P{6}; end;
if (nargs<7 | isempty(P{7})) Detrend = 'linear'; else Detrend = P{7}; end;
if (nargs<8 | isempty(P{8})) nTapers = 2*NW -1; else nTapers = P{8}; end;
if (nargs<9 | isempty(P{9})) FreqRange = [0 Fs/2]; else FreqRange = P{9}; end
if (nargs<10 | isempty(P{10})) TriggTimes = []; else TriggTimes = P{10}; TriggTimes=TriggTimes(:); end %<--- added parameter
if (nargs<11 | isempty(P{11})) ChannelPairs = []; else ChannelPairs = P{11};  end %<--- added parameter

% Now do some compuatations that are common to all spectrogram functions
if size(x,1)<size(x,2)
    x = x';
end
nChannels = size(x, 2);
nSamples = size(x,1);

if length(nOverlap)==1
    winstep = WinLength - nOverlap;
    % calculate number of FFTChunks per channel
    %remChunk = rem(nSamples-Window)
    nFFTChunks = max(1,round(((nSamples-WinLength)/winstep))); %+1  - is it ? but then get some error in the chunking in mtcsd... let's figure it later
    t = winstep*(0:(nFFTChunks-1))'/Fs;
else
    winstep = 0;
    nOverlap = nOverlap(nOverlap>WinLength/2 & nOverlap<nSamples-WinLength/2);
    nFFTChunks = length(nOverlap);
    t = nOverlap(:)/Fs; 
end 
%here is how welch.m of matlab does it:
% LminusOverlap = L-noverlap;
% xStart = 1:LminusOverlap:k*LminusOverlap;
% xEnd   = xStart+L-1;
% welch is doing k = fix((M-noverlap)./(L-noverlap)); why?
% turn this into time, using the sample frequency


% set up f and t arrays
if isreal(x)%~any(any(imag(x)))    % x purely real
	if rem(nFFT,2),    % nfft odd
		select = [1:(nFFT+1)/2];
	else
		select = [1:nFFT/2+1];
	end
	nFreqBins = length(select);
else
	select = 1:nFFT;
end
f = (select - 1)'*Fs/nFFT;
nFreqRanges = size(FreqRange,1);
%if (FreqRange(end)<Fs/2)
    if nFreqRanges==1
        select = find(f>FreqRange(1) & f<FreqRange(end));
        f = f(select);
        nFreqBins = length(select);
    else
        select=[];
        for i=1:nFreqRanges
            select=cat(1,select,find(f>FreqRange(i,1) & f<FreqRange(i,2)));
        end
        f = f(select);
        nFreqBins = length(select);
    end
%end

%-------------------------------------------------------------------------------%

function [M,IND] = combn(V,N)
% COMBN - all combinations of elements
%   M = COMBN(V,N) returns all combinations of N elements of the elements in
%   vector V. M has the size (length(V).^N)-by-N.
%
%   [M,I] = COMBN(V,N) also returns the index matrix I so that M = V(I).
%
%   V can be an array of numbers, cells or strings.
%
%   Example:
%     M = COMBN([0 1],3) returns the 8-by-3 matrix:
%       0     0     0
%       0     0     1
%       0     1     0
%       0     1     1
%       ...
%       1     1     1
%
%   All elements in V are regarded as unique, so M = COMBN([2 2],3) returns 
%   a 8-by-3 matrix with all elements equal to 2.
%
%   NB Matrix sizes increases exponentially at rate (n^N)*N.
% 
%   See also PERMS, NCHOOSEK
%        and ALLCOMB and PERMPOS on the File Exchange

% tested in Matlab R13, R14, 2010b
% version 4.2 (apr 2011)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History
% 1.1 updated help text
% 2.0 new faster algorithm
% 3.0 (aug 2006) implemented very fast algorithm
% 3.1 (may 2007) Improved algorithm Roger Stafford pointed out that for some values, the floor
% operation on floating points, according to the IEEE 754 standard, could return 
% erroneous values. His excellent solution was to add (1/2) to the values
% of A. 
% 3.2 (may 2007) changed help and error messages slightly
% 4.0 (may 2008) again a faster implementation, based on ALLCOMB, suggested on the
%     newsgroup comp.soft-sys.matlab on May 7th 2008 by "Helper". It was
%     pointed out that COMBN(V,N) equals ALLCOMB(V,V,V...) (V repeated N
%     times), ALLCMOB being faster. Actually version 4 is an improvement
%     over version 1 ... 
% 4.1 (jan 2010) removed call to FLIPLR, using refered indexing N:-1:1
%     (is faster, suggestion of Jan Simon, jan 2010), removed REPMAT, and
%     let NDGRID handle this
% 4.2 (apr 2011) corrrectly return a column vector for N = 1 (error pointed
%      out by Wilson).

error(nargchk(2,2,nargin)) ;

if isempty(V) || N == 0,
    M = [] ;
    IND = [] ;
elseif fix(N) ~= N || N < 1 || numel(N) ~= 1 ;
    error('combn:negativeN','Second argument should be a positive integer') ;
elseif N==1,
    % return column vectors
    M = V(:) ; 
    IND = (1:numel(V)).' ;
else
    % speed depends on the number of output arguments
    if nargout<2,
        M = local_allcomb(V,N) ;
    else
        % indices requested
        IND = local_allcomb(1:numel(V),N) ;
        M = V(IND) ;
    end
end

function Y = local_allcomb(X,N)
% See ALLCOMB, available on the File Exchange
if N>1
    % create a list of all possible combinations of N elements
    [Y{N:-1:1}] = ndgrid(X) ;
    % concatenate into one matrix, reshape into 2D and flip columns
    Y = reshape(cat(N+1,Y{:}),[],N) ;
else
    % no combinations have to be made
    Y = X(:) ;
end




