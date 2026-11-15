macro "Automatic Lipid Droplet Circles" {

    if (nImages == 0)
        exit("Open a fluorescence image first.");

    original   = getTitle();
    originalID = getImageID();

    // -------------------------------
    // USER SETTINGS
    // -------------------------------

    Dialog.create("Lipid Droplet Detection");

    Dialog.addNumber("Typical droplet diameter (pixels)", 20);
    Dialog.addNumber("Detection threshold (% of max, lower = more droplets)", 20);
    Dialog.addNumber("Minimum droplet area (pixels)", 40);
    Dialog.addNumber("Maximum droplet area (pixels)", 1200);
    Dialog.addNumber("Minimum circularity (0-1)", 0.55);
    Dialog.addCheckbox("Apply watershed", true);
    Dialog.addCheckbox("Draw ROIs as fitted circles", true);

    Dialog.show();

    diameter  = Dialog.getNumber();
    thrPct    = Dialog.getNumber();
    minArea   = Dialog.getNumber();
    maxArea   = Dialog.getNumber();
    minCirc   = Dialog.getNumber();
    watershed = Dialog.getCheckbox();
    fitCircle = Dialog.getCheckbox();

    if (diameter <= 0 || thrPct <= 0 || thrPct >= 100 ||
        minArea < 0 || maxArea <= minArea ||
        minCirc < 0 || minCirc > 1)
        exit("Check your parameters.");

    sigmaSmall = diameter / 4;   // fills ring-shaped droplets into blobs
    sigmaBig   = diameter;       // estimates the cytoplasm background

    // -------------------------------
    // SEGMENTATION (band-pass to keep droplet-sized objects only)
    // -------------------------------

    selectImage(originalID);

    // Duplicate ONLY the currently displayed channel / z / t plane
    if (is("hyperstack")) {
        Stack.getPosition(c, z, t);
        run("Duplicate...", "title=LD_MASK duplicate channels=" + c +
            " slices=" + z + " frames=" + t);
    } else {
        run("Duplicate...", "title=LD_MASK");
    }

    selectWindow("LD_MASK");
    run("32-bit");
    run("Duplicate...", "title=LD_BIG");

    selectWindow("LD_MASK");
    run("Gaussian Blur...", "sigma=" + sigmaSmall);

    selectWindow("LD_BIG");
    run("Gaussian Blur...", "sigma=" + sigmaBig);

    // Difference of Gaussians: droplets stay, smooth cytoplasm cancels out
    imageCalculator("Subtract 32-bit", "LD_MASK", "LD_BIG");

    selectWindow("LD_BIG");
    close();

    selectWindow("LD_MASK");
    getStatistics(area, mean, min, max);

    if (max <= 0) {
        close();
        exit("No droplet-sized structures found. Check droplet diameter.");
    }

    // Scale positive signal to 0-255, then threshold as % of max
    setMinAndMax(0, max);
    run("8-bit");

    setOption("BlackBackground", true);
    setThreshold(round(255 * thrPct / 100), 255);
    run("Convert to Mask");

    run("Fill Holes");

    if (watershed)
        run("Watershed");

    // -------------------------------
    // DETECT DROPLETS
    // -------------------------------

    roiManager("Reset");

    run("Analyze Particles...",
        "size=" + minArea + "-" + maxArea + " pixel" +
        " circularity=" + minCirc + "-1.00" +
        " show=Nothing add");

    n = roiManager("Count");

    selectWindow("LD_MASK");
    close();

    if (n == 0)
        exit("No droplets detected. Lower the threshold %, minimum circularity or minimum area.");

    // -------------------------------
    // MAKE ROIs ROUND AND NAME THEM
    // -------------------------------

    selectImage(originalID);

    for (i = 0; i < n; i++) {
        roiManager("Select", i);
        if (fitCircle) {
            run("Fit Circle");
            roiManager("Update");
        }
        roiManager("Rename", "LD_" + (i + 1));
    }

    // -------------------------------
    // MEASURE ORIGINAL INTENSITY
    // -------------------------------

    selectImage(originalID);

    run("Clear Results");
    run("Set Measurements...",
        "area mean min integrated centroid perimeter redirect=None decimal=3");

    roiManager("Deselect");
    roiManager("Measure");

    // -------------------------------
    // DISPLAY FINAL ROIs
    // -------------------------------

    roiManager("Deselect");
    roiManager("Show All without labels");
    roiManager("Set Color", "yellow");
    roiManager("Set Line Width", 1);

    print("Lipid droplet analysis completed.");
    print("Detected droplets: " + n);

    showMessage("Analysis completed",
        "Droplets detected: " + n +
        "\nMean intensity is in the Results table.");
}
