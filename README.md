# EC-SOFI
this matlab code listed here is for the publication 'Electrochemically controlled switching of dyes for enhanced super-resolution optical fluctuation imaging (EC-SOFI)' by Ying Yang, Yuanqing Ma, Richard D. Tilley, J. Justin Gooding. 
The code inclue:
1, various single molecule photophyisics analysis of the flurophore blinking of DNA origami as in figure 1, which uses Zeiss zen fitted STORM table.
2, SOFI image calculation from the raw zeiss file, where the time series were slitted into subsequences and SOFI calculation was performed for each subsequences and summed. The SOFI  calculation itself was conducted by following the bSOFI code published in (GitHub - lob-epfl/sofitool) with minor modificaations, which can be found in the published manuscript. 
3, Fourier ring correlation analysis for the SOFI of various order and STORM images. 
The raw data of the publication will be depoist into public data domain in DRYAD (https://datadryad.org/) shortly. 

