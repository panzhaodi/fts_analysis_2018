fts3g.py is the code for generating every single spectra. I ran the code three times last year for centered pixels, neighbor pixels, and second neighbors. I stored the three groups in three folders.

pyfts_package.py is the low-level fts analysis modules that will be referred by fts3g.py.

plot_corrected_bands.py   This code plots the bands corrected for the Lyot stop, Rayleigh-Jeans source, AOmega of the detector, and the nu^2 from the Mylar filter. All corrections are encapculated in fts_correction.txt. Brad generated this file (ask him if you have questions). So the corrected data should only be the spectrum of the detector (plus the lenses, lenslets, and metal mesh filters). 

plot_averaged_normalized.py As the name suggests, this code plot the averaged normalized spectra per wafer. It generates one plot per wafer and one curve for the centered pixels, neighboring pixels, and second neighbors separately.

plot_raw.py plots the raw spectra, unnormalized. It generates one plot per wafer, one curve for each detector. The 90s, 150s, and 220s are in different colors.

plot_raw_normalized.py is similar to plot_raw.py, but with normalized detector curves.

print_numbers.py will print the band stat numbers, like the band center, band edges and band widths. Note that this code is using a different correction function (nu^4), which only includes the nu^2 from the Mylar and the nu^2 from the blackbody. The Lyot stop and the AOmega of the detector are not corrected, because these two things are also there when the telescope observes the sky. So this is the effective bandwidth of the detector through the receiver (not only the detector). 

pixel_positions.pkl has all the pixel positions on a wafer for determining the neighbors and nearest neighbors.

fts_correction.txt has the correction factor for the spectra, and will be used for plot_corrected_bands.py. Its plot is correction_factor.png. 
