import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT

"""
Script to plot out (nicely) the PDF information for the different energy bins
for the directionality and multisite.

The plots generated are: 
    1. Time Residual comparison for Tl208 and B8 (multisite PDFs)
    2. Directionality cos(theta_gamma) vs b8 tRes PDF
    3. Energy distributions of PDF'd events
    4. Position distributions for entire dataset
        a. rho2-Z
        b. X-Y
    5. Radius vs Energy histogram
"""

def plot_pdf_analytics(pdf_energy_string):

    # use the string to define energy cuts
    if pdf_energy_string == "2.5_5.0":
        energy_low  = 2.5
        energy_high = 5.0
    if pdf_energy_string == "2.5_3.125":
        energy_low  = 2.5
        energy_high = 3.125
    if pdf_energy_string == "3.125_3.75":
        energy_low  = 3.125
        energy_high = 3.75
    if pdf_energy_string == "3.75_4.375":
        energy_low  = 3.75
        energy_high = 4.375
    if pdf_energy_string == "4.375_5.0":
        energy_low  = 4.375
        energy_high = 5.0
    
    # load the information
    pdf_B8_file    = ROOT.TFile.Open("../run_by_run_pdf/B8_Solar_Nue/total.root")
    pdf_Tl208_file = ROOT.TFile.Open("../run_by_run_pdf/Tl208/total.root")

    # extract the desired PDF histograms for re-plotting with matplotlib
    multisite_pdf_B8    = pdf_B8_file.Get(f"multi_{pdf_energy_string}")
    dir_pdf_B8          = pdf_B8_file.Get(f"directionality_{pdf_energy_string}")
    multisite_pdf_Tl208 = pdf_Tl208_file.Get(f"multi_{pdf_energy_string}")
    dir_pdf_Tl208       = pdf_Tl208_file.Get(f"directionality_{pdf_energy_string}")
    
    # get the bin edges and contents for each PDF
    multi_B8_bins      = [multisite_pdf_B8.GetBinLowEdge(i) for i in range(1, multisite_pdf_B8.GetNbinsX() + 2)]
    multi_B8_counts    = [multisite_pdf_B8.GetBinContent(i) / multisite_pdf_B8.Integral() for i in range(1, multisite_pdf_B8.GetNbinsX()+1)]
    multi_Tl208_bins   = [multisite_pdf_Tl208.GetBinLowEdge(i) for i in range(1, multisite_pdf_Tl208.GetNbinsX() + 2)]
    multi_Tl208_counts = [multisite_pdf_Tl208.GetBinContent(i) / multisite_pdf_Tl208.Integral() for i in range(1, multisite_pdf_Tl208.GetNbinsX()+1)]

    # create numpy arrays to save the directionality 2D PDF bins
    directionality_pdf = np.zeros((dir_pdf_B8.GetNbinsX(), dir_pdf_B8.GetNbinsY()))
    for i in range(1, dir_pdf_B8.GetNbinsX() + 1):
        for j in range(1, dir_pdf_B8.GetNbinsY() + 1):
            directionality_pdf[i-1, j-1] = dir_pdf_B8.GetBinContent(i,j)
    x_min = dir_pdf_B8.GetXaxis().GetXmin()
    x_max = dir_pdf_B8.GetXaxis().GetXmax()
    y_min = dir_pdf_B8.GetYaxis().GetXmin()
    y_max = dir_pdf_B8.GetYaxis().GetXmax()

    # get the event level info for plotting the position and energy distributions
    event_info_B8    = pdf_B8_file.Get("PDF_event_by_event_Information")
    event_info_Tl208 = pdf_Tl208_file.Get("PDF_event_by_event_Information")
    
    num_events_B8    = event_info_B8.GetEntries()
    num_events_Tl208 = event_info_Tl208.GetEntries()
    energy_B8        = np.zeros(num_events_B8)
    rho2_B8          = np.zeros(num_events_B8)
    z_B8             = np.zeros(num_events_B8)
    x_B8             = np.zeros(num_events_B8)
    y_B8             = np.zeros(num_events_B8)
    r_B8             = np.zeros(num_events_B8)
    energy_Tl208     = np.zeros(num_events_Tl208)
    rho2_Tl208       = np.zeros(num_events_Tl208)
    z_Tl208          = np.zeros(num_events_Tl208)
    x_Tl208          = np.zeros(num_events_Tl208)
    y_Tl208          = np.zeros(num_events_Tl208)
    r_Tl208          = np.zeros(num_events_Tl208)
    counter = 0
    for ientry in event_info_B8:
        e = ientry.energy
        if e >= energy_low and e < energy_high:
            energy_B8[counter] = e

            X = ientry.x
            Y = ientry.y
            Z = ientry.z
            rho2_B8[counter] = (X**2 + Y**2) / 6000**2
            z_B8[counter]    = Z
            x_B8[counter]    = X
            y_B8[counter]    = Y
            r_B8[counter]    = np.sqrt(X**2 + Y**2 + Z**2)
        counter += 1
    counter = 0
    for ientry in event_info_Tl208:
        e = ientry.energy
        if e >= energy_low and e < energy_high:
            energy_Tl208[counter] = e

            X = ientry.x
            Y = ientry.y
            Z = ientry.z
            rho2_Tl208[counter] = (X**2 + Y**2) / 6000**2
            z_Tl208[counter]    = Z
            x_Tl208[counter]    = X
            y_Tl208[counter]    = Y
            r_Tl208[counter]    = np.sqrt(X**2 + Y**2 + Z**2)
        counter += 1
    
    # remove all the zeros (if any) left in the arrays
    rho2_B8 = rho2_B8[np.nonzero(rho2_B8)]
    z_B8 = z_B8[np.nonzero(z_B8)]
    x_B8 = x_B8[np.nonzero(x_B8)]
    y_B8 = y_B8[np.nonzero(y_B8)]
    energy_B8 = energy_B8[np.nonzero(energy_B8)]
    r_B8 = r_B8[np.nonzero(r_B8)]

    rho2_Tl208 = rho2_Tl208[np.nonzero(rho2_Tl208)]
    z_Tl208 = z_Tl208[np.nonzero(z_Tl208)]
    x_Tl208 = x_Tl208[np.nonzero(x_Tl208)]
    y_Tl208 = y_Tl208[np.nonzero(y_Tl208)]
    energy_Tl208 = energy_Tl208[np.nonzero(energy_Tl208)]
    r_Tl208 = r_Tl208[np.nonzero(r_Tl208)]

    fig, axes = plt.subplots(nrows = 3, ncols = 3, figsize = (25, 20))

    # multisite PDF
    axes[0,0].step(multi_B8_bins[:-1], multi_B8_counts,  where = 'post', label = r"$^8B$")
    axes[0,0].step(multi_B8_bins[:-1], multi_Tl208_counts,  where = 'post', label = r"$^{208}Tl$")
    axes[0,0].legend(frameon=False, fontsize = 20)
    axes[0,0].set_xlabel("Time Residual (ns)", fontsize = 20)
    axes[0,0].set_ylabel("Normalised Counts per 1 ns Bin", fontsize = 20)
    axes[0,0].set_xlim((-5, 100))

    # directionality PDF
    axes[0,1].imshow(directionality_pdf.T, origin = 'lower', aspect = 'auto', extent = [x_min, x_max, y_min, y_max])
    axes[0,1].set_xlabel("Time Residual (ns)", fontsize = 20)
    axes[0,1].set_ylabel(r"$cos(\theta _\gamma)$", fontsize = 20)
    axes[0,1].set_xlim((-10, 80))

    # energy distribution
    width = 0.05
    energy_bins = np.arange(2.5, 5.0 + width, width)
    axes[0,2].hist(energy_B8, bins = energy_bins, density = True, histtype = "step", label = r"$^8B$")
    axes[0,2].hist(energy_Tl208, bins = energy_bins, density = True, histtype = "step", label = r"$^{208}Tl$")
    axes[0,2].set_xlabel("Reconstructed Energy (MeV)", fontsize = 20)
    axes[0,2].set_ylabel(f"Normalised Counter per {width} MeV", fontsize = 20)
    axes[0,2].legend(frameon=False, fontsize = 20)
    
    # energy vs radius
    r_bins = np.linspace(0, 6000, 50)
    axes[1,0].hist2d(r_B8, energy_B8, bins = [r_bins, energy_bins], cmin = 1e-4, density = True)
    axes[1,0].set_xlabel("Reconstructed Radius (mm)", fontsize = 20)
    axes[1,0].set_ylabel("Reconstructed Energy (MeV)", fontsize = 20)
    axes[1,0].plot([], [], linestyle = "", label = r"$^8B$")
    axes[1,0].legend(frameon=False, fontsize = 20)

    axes[1,1].hist2d(r_Tl208, energy_Tl208, bins = [r_bins, energy_bins], cmin = 1e-4, density = True)
    axes[1,1].set_xlabel("Reconstructed Radius (mm)", fontsize = 20)
    axes[1,1].set_ylabel("Reconstructed Energy (MeV)", fontsize = 20)
    axes[1,1].plot([], [], linestyle = "", label = r"$^{208}Tl$")
    axes[1,1].legend(frameon = False, fontsize = 20)

    # rho2_z position distribution
    z_bins    = np.linspace(-6000, 6000, 100)
    rho2_bins = np.linspace(0, 1, 100)
    axes[1,2].hist2d(rho2_B8, z_B8, bins = [rho2_bins, z_bins], density = True, cmin = 1e-4)
    axes[1,2].plot([], [], linestyle = "", label = r"$^8B$")
    axes[1,2].set_xlabel(r"$\left(\frac{\rho}{\rho_{AV}}\right)^2$", fontsize = 20)
    axes[1,2].set_ylabel("Reconstructed Z (mm)", fontsize = 20)
    axes[1,2].legend(frameon=False, fontsize = 20)

    axes[2,0].hist2d(rho2_Tl208, z_Tl208, bins = [rho2_bins, z_bins], density = True, cmin = 1e-4)
    axes[2,0].plot([], [], linestyle = "", label = r"$^{208}Tl$")
    axes[2,0].set_xlabel(r"$\left(\frac{\rho}{\rho_{AV}}\right)^2$", fontsize = 20)
    axes[2,0].set_ylabel("Reconstructed Z (mm)", fontsize = 20)
    axes[2,0].legend(frameon=False, fontsize = 20)

    # X-Y position plot
    axes[2,1].hist2d(x_B8, y_B8, bins = [z_bins, z_bins], density = True, cmin = 1e-10)
    axes[2,1].set_xlabel("Reconstructed X (mm)", fontsize = 20)
    axes[2,1].set_ylabel("Reconstructed Y (mm)", fontsize = 20)
    axes[2,1].plot([], [], linestyle = "", label = r"$^8B$")
    axes[2,1].legend(frameon = False, fontsize = 20)

    axes[2,2].hist2d(x_Tl208, y_Tl208, bins = [z_bins, z_bins], density = True, cmin = 1e-10)
    axes[2,2].set_xlabel("Reconstructed X (mm)", fontsize = 20)
    axes[2,2].set_ylabel("Reconstructed Y (mm)", fontsize = 20)
    axes[2,2].plot([], [], linestyle = "", label = r"$^{208}Tl$")
    axes[2,2].legend(frameon = False, fontsize = 20)
    
    fig.tight_layout()
    plt.subplots_adjust(top=0.95)
    plt.suptitle(f"PDF Information: {pdf_energy_string} MeV", fontsize = 30)
    plt.savefig(f"../plots/pdf_analytics/{pdf_energy_string}.pdf")
plot_pdf_analytics("2.5_5.0")
plot_pdf_analytics("2.5_3.125")
plot_pdf_analytics("3.125_3.75")
plot_pdf_analytics("3.75_4.375")
plot_pdf_analytics("4.375_5.0")