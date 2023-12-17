function dataOut = receiveSignal(txSignal, modOrder, M, numChannels, powersVec, sgm, U, varNoise, typeDetector, uRot)
    useRot = ~isnan(uRot);

    rxSignal = U' * txSignal;

    for k = 1:numChannels
        if(strcmp(typeDetector, 'MMSE'))
            outputSymbols(k, 1) = rxSignal(k) * (sgm(k, k) * powersVec(k)) / (sgm(k, k)^2 * powersVec(k)^2 + varNoise);
        elseif(strcmp(typeDetector, 'ZF'))
            outputSymbols(k, 1) = rxSignal(k) / (sgm(k, k) * powersVec(k));
        end
    end

    if useRot   
        outputSymbols = uRot' * outputSymbols;
    else
        outputSymbols = outputSymbols;
    end
        
    outputData = qamdemod(outputSymbols, modOrder);

    dataOut = de2bi(outputData, M);
end

