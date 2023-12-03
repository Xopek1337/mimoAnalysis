clc; clear;
%close all

txAntennasNum = 2;
rxAntennasNum = 4;

numChannels = 1;

M = 4;
modOrder = 2^M;

typeDetector = 'ZF';

SNR_dB = (0:1:30);

varNoise = 1;

nRealiz = 10000;

nErr = zeros(length(SNR_dB), nRealiz);

for i = 1:txAntennasNum
    for j = 1:rxAntennasNum
        uRot(i, j) = (1/sqrt(numChannels)) * exp( ((1i*2*pi)/numChannels) * (j-1) * (i-1) );
    end
end

h = (randn(txAntennasNum, rxAntennasNum) + 1i*randn(txAntennasNum, rxAntennasNum)) *(1/sqrt(2));

[U, sgm, V] = svd(h);

sgm = sgm(1:numChannels, 1:numChannels);

U = U(:, 1:numChannels);
V = V(:, 1:numChannels);

for i = 1:length(SNR_dB)
    for j = 1:nRealiz
        SNR = 10^(SNR_dB(i)/10);

        Pin = SNR * varNoise;

        powersVec = sqrt(Pin/numChannels)*ones(1, numChannels);
        powersMat = diag(powersVec);

        inputData = randi([0 1], numChannels, M);

        rxSignal = signalTransmit(inputData, modOrder, numChannels, h, powersMat, U, V, SNR_dB(i), Pin, varNoise, NaN);
        
        for k = 1:numChannels
            if(strcmp(typeDetector, 'MMSE'))
                outputSymbols(k, 1) = rxSignal(k) / (sgm(k, k) * powersVec(k) + varNoise);
            elseif(strcmp(typeDetector, 'ZF'))
                outputSymbols(k, 1) = rxSignal(k) / (sgm(k, k) * powersVec(k));
            end
        end

        outputSymbols2 = outputSymbols / sqrt(Pin);
        
        outputData = qamdemod(outputSymbols2, modOrder, 'UnitAveragePower', true);

        dataOut = de2bi(outputData, M);
    
        [nErrors, ~] = biterr(inputData, dataOut);

        nErr(i, j) = nErrors;
    end
end

for i = 1:length(SNR_dB)
    for j = 1:nRealiz
        SNR = 10^(SNR_dB(i)/10);

        Pin = SNR * varNoise;

        powersVec = sqrt(Pin/numChannels)*ones(1, numChannels);
        powersMat = diag(powersVec);

        inputData = randi([0 1], numChannels, M);

        rxSignal = signalTransmit(inputData, modOrder, numChannels, h, powersMat, U, V, SNR_dB(i), Pin, varNoise, uRot);
        
        for k = 1:numChannels
            if(strcmp(typeDetector, 'MMSE'))
                outputSymbols(k, 1) = rxSignal(k) / (sgm(k, k) * powersVec(k) + varNoise);
            elseif(strcmp(typeDetector, 'ZF'))
                outputSymbols(k, 1) = rxSignal(k) / (sgm(k, k) * powersVec(k));
            end
        end

        outputSymbols = uRot(1:numChannels, 1:numChannels)' * outputSymbols;

        outputSymbols = outputSymbols / sqrt(Pin);
        
        outputData = qamdemod(outputSymbols, modOrder, 'UnitAveragePower' , true);

        dataOut = de2bi(outputData, M);
    
        [nErrors, ~] = biterr(inputData, dataOut);

        nErrWRotate(i, j) = nErrors;
    end
end

ber = sum(nErr,2)./(nRealiz*numChannels*M);

ber2 = sum(nErrWRotate,2)./(nRealiz*numChannels*M);

figure;
semilogy(SNR_dB, ber);
hold on;
semilogy(SNR_dB, ber2); grid on;
hold off;
title('16-QAM');
legend('Без матрицы поворота', 'С матрицей поворота');