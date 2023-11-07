function rxSignal = formSignal(inputData, modOrder, numSymbs, h, powersMat, U, V, snr, uRot)
    useRot = ~isnan(uRot);
    
    dataSym = bi2de(inputData);

    inputSymbols = qammod(dataSym, modOrder, 'UnitAveragePower' , true);

    if useRot    
        txSignal = h * powersMat *  V *  uRot * inputSymbols;
    else
        txSignal = h * V * powersMat * inputSymbols;
    end
        
    varSignal = var(txSignal);
    varNoise = varSignal*10^(-snr/10);
    txSignal = txSignal + sqrt(varNoise/2)*(randn(size(txSignal))+1i*randn(size(txSignal)));
          
    rxSignal = U' * txSignal;

end

