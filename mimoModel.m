clc; clear;
%close all

txAntennasNum = 8;
rxAntennasNum = txAntennasNum;

numSymbs = txAntennasNum;

M = 4;
modOrder = 2^M;

snr = (-10:1:45);

Pin = 1;

powersVec = sqrt(Pin/numSymbs)*ones(1, numSymbs);

powersMat = diag(powersVec);

nRealiz = 15000;

nErr = zeros(length(snr), nRealiz);

for i = 1:txAntennasNum
    for j = 1:rxAntennasNum
        uRot(i, j) = (1/sqrt(numSymbs)) * exp( ((1i*2*pi)/numSymbs) * (j-1) * (i-1) );
    end
end

h = reshape((randn(txAntennasNum * rxAntennasNum, 1) + 1i*randn(txAntennasNum * rxAntennasNum, 1))*(1/sqrt(2)), [txAntennasNum, rxAntennasNum]);

[U, sgm, V] = svd(h);

inputData = randi([0 1], numSymbs, M);

for i = 1:length(snr)
    for j = 1:nRealiz
        rxSignal = formSignal(inputData, modOrder, numSymbs, h, powersMat, U, V, snr(i), NaN);
        
        for k = 1:numSymbs
            outputSymbols(k, 1) = rxSignal(k) / sgm(k, k) / powersVec(k);
        end
        
        outputData = qamdemod(outputSymbols, modOrder, 'UnitAveragePower' , true);

        dataOut = de2bi(outputData, M);
    
        [nErrors, ~] = biterr(inputData, dataOut);

        nErr(i, j) = nErrors;
    end
end

for i = 1:length(snr)
    for j = 1:nRealiz
        rxSignal = formSignal(inputData, modOrder, numSymbs, h, powersMat, U, V, snr(i), uRot);
        
        for k = 1:numSymbs
            outputSymbols(k, 1) = rxSignal(k) / sgm(k, k) / powersVec(k);
        end

        outputSymbols = uRot' * outputSymbols;
        
        outputData = qamdemod(outputSymbols, modOrder, 'UnitAveragePower' , true);

        dataOut = de2bi(outputData, M);
    
        [nErrors, ~] = biterr(inputData, dataOut);

        nErrWRotate(i, j) = nErrors;
    end
end

ber = sum(nErr,2)./(nRealiz*numSymbs*M);

ber2 = sum(nErrWRotate,2)./(nRealiz*numSymbs*M);

semilogy(snr, ber);
hold on;
semilogy(snr, ber2); grid on;
hold off;
title('16-QAM');
legend('Без матрицы поворота', 'С матрицей поворота');
