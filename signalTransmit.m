function rxSignal = signalTransmit(inputData, modOrder, numChannels, h, powersMat, U, V, Pin, varNoise, uRot)
    useRot = ~isnan(uRot);
    
    dataSym = bi2de(inputData);

    inputSymbols = qammod(dataSym, modOrder, 'UnitAveragePower', true) * sqrt(Pin);

    if useRot    
        txSignal = h *  V * powersMat * uRot(1:numChannels, 1:numChannels) * inputSymbols;
    else
        txSignal = h * V * powersMat * inputSymbols;
    end

    noise =  sqrt(varNoise)*(randn(size(txSignal))+1i*randn(size(txSignal)));
    txSignal = txSignal + noise;
          
    rxSignal = U' * txSignal;
end

