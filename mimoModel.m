clc; clear;
%close all

txAntennasNum = 4;
rxAntennasNum = 4;

numChannels = 2;

M = 4;
modOrder = 2^M;

typeDetector = 'ZF';

seed = 1;

SNR_dB = (-10:1:20);

varNoise = 1;

nRealiz = 30000;

nErr = zeros(length(SNR_dB), nRealiz);
nErrWRotate = nErr;

for i = 1:numChannels
    for j = 1:numChannels
        uRot(i, j) = (1/sqrt(numChannels)) * exp( ((1i*2*pi)/numChannels) * (i-1) * (j-1) );
    end
end

randn('state', seed);

h = zeros(txAntennasNum, rxAntennasNum);

while(abs(mean(mean(abs(h.^2))) - 1) > 0.001)
    h = (randn(txAntennasNum, rxAntennasNum) + 1i*randn(txAntennasNum, rxAntennasNum)) *(1/sqrt(2));
end

[U, sgm, V] = svd(h);

sgm = sgm(1:numChannels, 1:numChannels);

U = U(:, 1:numChannels);
V = V(:, 1:numChannels);

percArray   = 0;
disp([num2str(percArray(end)) ' %']);

for i = 1:length(SNR_dB)
    SNR = 10^(SNR_dB(i)/10);

    Pin = SNR * varNoise;

    powersVec = sqrt(Pin/numChannels)*ones(1, numChannels);
    powersMat = diag(powersVec);

    for j = 1:nRealiz
        % Показывает прогресс выполнения программы
        indPrt = round(((j+(i-1)*nRealiz)/(nRealiz*length(SNR_dB))) * 1e2);
        perc   = indPrt; 
        if perc ~= percArray(end) 
            percArray(end + 1) = perc;
            disp([num2str(percArray(end)) ' %']);
        end

        inputData = randi([0 1], numChannels, M);

        inputSymbols = formSymbols(inputData, modOrder);

        txSignal = transmitSignal(inputSymbols, h, powersMat, V, varNoise, NaN);
        dataOut = receiveSignal(txSignal, modOrder, M, numChannels, powersVec, sgm, U, varNoise, typeDetector, NaN);

        txSignalRot = transmitSignal(inputSymbols, h, powersMat, V, varNoise, uRot);
        dataOutRot = receiveSignal(txSignalRot, modOrder, M, numChannels, powersVec, sgm, U, varNoise, typeDetector, uRot);

        [nErrors, ~] = biterr(inputData, dataOut);
        [nErrorsWRotate, ~] = biterr(inputData, dataOutRot);

        nErr(i, j) = nErrors;
        nErrWRotate(i, j) = nErrorsWRotate;
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