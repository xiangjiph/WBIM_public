yGalvoWaveform = hSI.hWaveformManager.scannerAO.ao_volts_raw.G(:,2);
yGalvoSampleRate = hSI.hWaveformManager.scannerAO.sampleRates.G;

rs = dabs.resources.ResourceStore();
galvo = rs.filterByName('Galvo');
galvo.optimizeWaveformIteratively(yGalvoWaveform,yGalvoSampleRate);