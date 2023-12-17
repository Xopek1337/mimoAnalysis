function inputSignal = formSymbols(inputData, modOrder, numChannels, h, powersMat, U, V, Pin, varNoise, uRot)
    dataSym = bi2de(inputData);

    inputSignal = qammod(dataSym, modOrder);
end

