function rxSignal = signalTransmit(inputData, modOrder, numSymbs, h, powersMat, U, V, Pin, varNoise, uRot)
    useRot = ~isnan(uRot);
    
    dataSym = bi2de(inputData);

    inputSymbols = qammod(dataSym, modOrder, 'UnitAveragePower', true) * sqrt(Pin);

    if useRot    
        txSignal = h *  V * powersMat * uRot(1:numSymbs, 1:numSymbs) * inputSymbols;
    else
        txSignal = h * V * powersMat * inputSymbols;
    end

    noise =  sqrt(varNoise)*(randn(size(txSignal))+1i*randn(size(txSignal)));
    txSignal = txSignal + noise;
          
    rxSignal = U' * txSignal;
end

