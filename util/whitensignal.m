%function [y, ARmodel] = WhitenSignal(x, window,CommonAR,ARmodel, ARorder, FilterFun)
% whitens the signal 
% if window specified will recompute the model in each window of that size
% (window is in samples,e,g, 300sec*1250 samples
% if CommonAR is set to 1, then will use model from first channel for all
% if ARmodel is specified - use it, not compute fromthe data
% output optionaly the ARmodel for use on the other data to be on the same scale
% new parameter: FilterFun : filter or filtfilt
% NB!!!!!!!!!!!!!!!
% For spectrum filter is better , as filtfilt filters by
% square of the AR model, need to find a way to fix that
% For phase filtfilt is definitely better. 
% for now - default is filter (to make spectrum reasonable)
function [y, A] = WhitenSignal(x,varargin)

[window,CommonAR, ARmodel,ArOrder, FilterFun] = DefaultArgs(varargin,{[],1,[],2,'filter'});
% warning('Please, check help to understand what filter settings are good for what purposes');

Trans = 0;
if size(x,1)<size(x,2)
    x = x';
    Transf =1;
end
[nT nCh]  = size(x);
y = zeros(nT,nCh);

if isempty(window)
    seg = [1 nT];
    nwin=1;
else
    nwin = floor(nT/window)+1;
    seg = repmat([1 window],nwin,1)+repmat([0:nwin-1]'*window,1,2);
    if nwin*window>nT
        seg(end,2) =nT;
    end   
end

for j=1:nwin
    if ~isempty(ARmodel) 
        A = ARmodel;
        for i=1:nCh
            y(seg(j,1):seg(j,2),i) = RemoveAR(x(seg(j,1):seg(j,2),i),  A, FilterFun);
          
        end
    else
        if CommonAR % meaning common model for all channels and segments!!! 
            for i=1:nCh
                if  j==1 & i==1
                    A = arburg(x(seg(j,1):seg(j,2),i),ArOrder);
                end
                y(seg(j,1):seg(j,2),i) = RemoveAR(x(seg(j,1):seg(j,2),i),  A, FilterFun);
                
            end
        else
            for i=1:nCh
                A = arburg(x(seg(j,1):seg(j,2),i),ArOrder);
                y(seg(j,1):seg(j,2),i) = RemoveAR(x(seg(j,1):seg(j,2),i),  A,  FilterFun);
            end
        end
    end
end

if Trans
    y =y';
end
return

function cs = RemoveAR(s, a,ffun);
%cs = filter( a, 1, s);
cs = feval(ffun, a,1,s);
return



%%%%%%%%%%%%%%%%%%%%
A=[1 -2.7607 3.8106 -2.6535 0.9238];
% AR(4) coefficients
e = 0.2*randn(1024,1);
y=filter(1,A,e);
% Filter a white noise input to create AR(4) process
ar_coeffs=arburg(y,4);
%compare the results in ar_coeffs to the vector A

wy = filter(ar_coeffs, 1, y);
wy1 = filtfilt(ar_coeffs, 1, y);
plot([e wy wy1]);


x1 = x(1:1250*50);
Aburg = arburg(x1,4);

wx1 = filter(Aburg,1, x1);
wx2 = filtfilt(Aburg, 1, x1);
wx3 = Filter0(Aburg,  x1);


[y1 f1 phi1] = mtchd([x1 wx1 wx2 wx3],2^11,1250,2^10,[],2,'constant',[],[1 50]);

[y2 f2 phi2] = mtchd([x1 wx1 wx2 wx3],2^8,1250,2^7,[],2,'constant',[],[50 300]);
pow1 = MatDiag(y1);pow2 = MatDiag(y2);

figure;
clf
subplot(211)
hold on
plot(f1,log10(sq(pow1))); 
plot(f2, log10(sq(pow2)));
legend('original','filter','filtfilt','Filter0');

subplot(212)
hold on
plot(f1,sq(phi1(:,1,2:4))*180/pi); 
plot(f2, sq(phi2(:,1,2:4))*180/pi);

legend('filter','filtfilt','Filter0');


