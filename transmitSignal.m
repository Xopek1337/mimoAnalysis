function txSignal = transmitSignal(inputSymbols, h, powersMat, V, varNoise, uRot)
    useRot = ~isnan(uRot);
    
    if useRot    
        txSignal = h *  V * powersMat * uRot * inputSymbols;
    else
        txSignal = h * V * powersMat * inputSymbols;
    end

    noise =  sqrt(varNoise/2)*(randn(size(txSignal))+1i*randn(size(txSignal)));
    txSignal = txSignal + noise;
end

