import ROOT
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def create_roc(likelihood_signal, likelihood_background):

	# bin the data 
	min_vbb = np.amin(likelihood_signal)
	max_vbb = np.amax(likelihood_signal)
	min_cobalt = np.amin(likelihood_background)
	max_cobalt = np.amax(likelihood_background)
	if min_vbb < min_cobalt:
		bin_min = min_vbb
	else: 
		bin_min = min_cobalt
	if max_vbb > max_cobalt:
		bin_max = max_vbb
	else: 
		bin_max = max_cobalt
	print("Max Bin: {}\nMin Bin: {}".format(bin_max, bin_min))
	steps = 100#int((bin_max-bin_min)/0.2)
	print(steps)
	binning = np.linspace(bin_min, bin_max, steps)
	sig_counts, sig_bins = np.histogram(likelihood_signal, density=True, bins = binning)
	back_counts, back_bins = np.histogram(likelihood_background, density=True, bins = binning)

	max_hist = np.argmax(sig_counts)
	print("Number entries in signal hist: {}\nNum entries in background his: {}".format(np.size(likelihood_signal), np.size(likelihood_background)))
	# plt.figure()
	# plt.hist(likelihood_signal, density=True, histtype="step", bins = binning, label = "signal", linewidth=2)
	# plt.hist(likelihood_background, density=True, histtype="step", bins = binning, label = "background", linewidth=2)
	# plt.title("Log-Likelihood Ratio")
	# plt.legend()
	
	signal_acceptance = []
	background_acceptance = []
	
    # make a cut at each bin 
	tot_sig = np.sum(sig_counts)
	tot_back = np.sum(back_counts)
	print("LEN BINNING: ", len(binning[:-1])) 
	for i in range(len(binning[:-1])):
		# calculate the % of bin counts to the LEFT of the cut (hopefully the same as integrating it!)
		signal_accepted = np.sum(sig_counts[i:] * np.diff(sig_bins)[0]) 
		background_accepted = np.sum(back_counts[i:] * np.diff(back_bins)[0]) 
		signal_acceptance.append(signal_accepted)
		background_acceptance.append(background_accepted)

	return signal_acceptance, background_acceptance

def analysis_output_plot():
    """
    Function creates an absolutely massive (5 x 4) plot array showing: 
        1. dlog(L) for each energy bin
        2. ROC for multisite for each energy bin
        3. cos(theta_sun) for each energy bin
        4. Scatter of (cos(theta_sun), dlog(L)) for each energy bin 
    """

    # create the mega grid of plots
    fig, axes = plt.subplots(nrows = 3, ncols = 5, figsize = (30, 20))
    font_s    = 30
    
    bins_dir = np.linspace(-1, 1, 40)
    # load each of the test sets
    working_dir = "/data/snoplus3/hunt-stokes/multisite_directionality_0p6gL_refactored"
    test_B8     = ROOT.TFile.Open(working_dir + "/run_by_run_test/B8_Solar_Nue/total.root")
    test_Tl208  = ROOT.TFile.Open(working_dir + "/run_by_run_test/Tl208/total.root")

    # get each of the histograms
    histo_names = ["2p5_5p0", "2p5_3p125", "3p125_3p75", "3p75_4p375", "4p375_5p0"]
    plot_titles = ["2.5 --> 5.0 MeV", "2.5 --> 3.125 MeV", "3.125 --> 3.75 MeV", "3.75 --> 4.375 MeV", "4.375 --> 5.0 MeV"]
    plot_col    = 0
    for name in histo_names:
        
        dlogL_B8      = []
        dlogL_Tl208   = []
        cos_sun_B8    = []
        cos_sun_Tl208 = []
        roc_multi     = []

        # fill in the plots for each histogram in turn
        tree_b8    = test_B8.Get(name)
        tree_tl208 = test_Tl208.Get(name) 
        for ientry in tree_b8:
            dlogL_B8.append(ientry.dlogL)
            cos_sun_B8.append(ientry.cos_theta_sun)
        for ientry in tree_tl208:
            dlogL_Tl208.append(ientry.dlogL)
            cos_sun_Tl208.append(ientry.cos_theta_sun)
        
        # compute the ROC curve
        signal_roc, background_roc = create_roc(dlogL_B8, dlogL_Tl208)
		
        # fill the plot!
        axes[0, plot_col].hist(dlogL_B8, bins = 100, density = True, histtype = "step", label = r"$^8B$")
        axes[0, plot_col].hist(dlogL_Tl208, bins = 100, density = True, histtype = "step", label = r"$^{208}Tl$")
        axes[0, plot_col].legend(frameon=False, fontsize = font_s, loc = "upper left")
        axes[0, plot_col].set_xlabel(r"$\Delta log(\mathcal{L})$", fontsize = font_s)
        axes[0, plot_col].set_ylabel("Normalised Counts", fontsize = font_s)
        axes[0, plot_col].set_title(plot_titles[plot_col], fontsize = font_s)
        
        axes[1, plot_col].plot(signal_roc, background_roc, linewidth = 2)
        axes[1, plot_col].set_xlabel("Signal Acceptance", fontsize = font_s)
        axes[1, plot_col].set_ylabel("Background Acceptance", fontsize = font_s)
        axes[1, plot_col].grid("both")

        axes[2, plot_col].hist(cos_sun_B8, bins = bins_dir, density = True, histtype = "step", label = r"$^8B$")
        axes[2, plot_col].hist(cos_sun_Tl208, bins = bins_dir, density = True, histtype = "step", label = r"$^{208}Tl$")
        axes[2, plot_col].legend(frameon=False, fontsize = font_s, loc = "upper left")
        axes[2, plot_col].set_xlabel(r"$cos(\theta_{sun})$", fontsize = font_s)
        axes[2, plot_col].set_ylabel("Normalised Counts", fontsize = font_s)

        # axes[3, plot_col].scatter(dlogL_B8, cos_sun_B8, alpha = 0.1, s = 1, color = "red", label = r"$^8B$")
        # axes[3, plot_col].scatter(dlogL_Tl208, cos_sun_Tl208, alpha = 0.1, s = 1, color = "blue", label = r"$^{208}Tl$")
        # axes[3, plot_col].legend(frameon = False, fontsize = font_s, loc = "upper left")

        plot_col += 1
    fig.tight_layout()
    plt.savefig(working_dir + "/plots/test.pdf")

analysis_output_plot()
