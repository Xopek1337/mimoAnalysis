clc; clear;
close all

txAntennasNum = 4;
rxAntennasNum = txAntennasNum;

numSymbs = txAntennasNum;

M = 2;
modOrder = 2^M;

snr = 100;

Pin = 1;

powersVec = sqrt(Pin/numSymbs)*ones(1, numSymbs);

powersMat = diag(powersVec);

inputData = randi([0 modOrder-1], numSymbs, 1); 

inputSymbols = qammod(inputData, modOrder, 'UnitAveragePower' , true);

for i = 1:txAntennasNum
    for k = 1:rxAntennasNum
        Urot(i,k) = (1/sqrt(numSymbs)) * exp( ((1i*2*pi)/numSymbs) * (k-1) * (i-1) );
    end
end

h = reshape((randn(txAntennasNum * rxAntennasNum, 1) + 1i*randn(txAntennasNum * rxAntennasNum, 1))*(1/sqrt(2)), [txAntennasNum, rxAntennasNum]);

[U, sgm, V] = svd(h);

txSignal = h * V * powersMat * inputSymbols;

varSignal = var(txSignal);
varNoise = varSignal*10^(-snr/10);
txSignal = txSignal + sqrt(varNoise/2)*(randn(size(txSignal))+1i*randn(size(txSignal)));

rxSignal = U' * txSignal;

for i = 1:numSymbs
    outputSymbols(i, 1) = rxSignal(i) / sgm(i,i) / powersVec(i);
end

outputData = qamdemod(outputSymbols, modOrder, 'UnitAveragePower' , true);

disp(isequal(inputData, outputData));