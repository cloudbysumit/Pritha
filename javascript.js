macro "Automatic Lipid Droplet Circles (ring based)" {

    if (nImages == 0)
        exit("Open a fluorescence image first.");

    originalID = getImageID();

    // -------------------------------
    // USER SETTINGS
    // -------------------------------

    Dialog.create("Lipid Droplet Detection (ring based)");

    Dialog.addNumber("Blur sigma (pixels)", 1);
    Dialog.addNumber("Local threshold radius (pixels, ~ droplet diameter)", 20);
    Dialog.addNumber("Local threshold offset (higher = stricter)", 0);
    Dialog.addNumber("Minimum droplet area (pixels)", 40);
    Dialog.addNumber("Maximum droplet area (pixels)", 1200);
    Dialog.addNumber("Minimum circularity (0-1)", 0.70);
    Dialog.addNumber("Min contrast (inside / surround)", 1.15);
    Dialog.addCheckbox("Apply watershed", true);
    Dialog.addCheckbox("Draw ROIs as fitted circles", false);

    Dialog.show();

    sigma     = Dialog.getNumber();
    locRadius = Dialog.getNumber();
    locOffset = Dialog.getNumber();
    minArea   = Dialog.getNumber();
    maxArea   = Dialog.getNumber();
    minCirc   = Dialog.getNumber();
    minContr  = Dialog.getNumber();
    watershed = Dialog.getCheckbox();
    fitCircle = Dialog.getCheckbox();

    if (locRadius <= 0 || minArea < 0 || maxArea <= minArea ||
        minCirc < 0 || minCirc > 1 || sigma < 0)
        exit("Check your parameters.");

    // -------------------------------
    // SEGMENTATION
    // -------------------------------

    selectImage(originalID);

    // Duplicate ONLY the currently displayed channel / z / t plane
    if (is("hyperstack")) {
        Stack.getPosition(c, z, t);
        run("Duplicate...", "title=LD_A duplicate channels=" + c +
            " slices=" + z + " frames=" + t);
    } else {
        run("Duplicate...", "title=LD_A");
    }

    selectWindow("LD_A");
    run("8-bit");
    if (sigma > 0)
        run("Gaussian Blur...", "sigma=" + sigma);

    // Copy B: local threshold (keeps rings that are brighter than their surroundings)
    run("Duplicate...", "title=LD_B");
    selectWindow("LD_B");
    run("Auto Local Threshold",
        "method=Mean radius=" + locRadius +
        " parameter_1=" + locOffset + " parameter_2=0 white");

    // Copy A: global floor, so dark background noise is ignored
    selectWindow("LD_A");
    setOption("BlackBackground", true);
    setAutoThreshold("Triangle dark");
    run("Convert to Mask");

    // Keep only pixels that pass BOTH tests
    imageCalculator("AND create", "LD_A", "LD_B");
    rename("LD_MASK");

    selectWindow("LD_A"); close();
    selectWindow("LD_B"); close();

    selectWindow("LD_MASK");
    setOption("BlackBackground", true);

    run("Close-");        // close small gaps in rings
    run("Fill Holes");    // closed rings -> solid discs
    run("Open");          // remove thin leftovers

    if (watershed)
        run("Watershed");

    // -------------------------------
    // DETECT ROUND OBJECTS
    // -------------------------------

    roiManager("Reset");

    run("Analyze Particles...",
        "size=" + minArea + "-" + maxArea + " pixel" +
        " circularity=" + minCirc + "-1.00" +
        " show=Nothing exclude add");

    selectWindow("LD_MASK");
    close();

    n = roiManager("Count");
    if (n == 0)
        exit("No droplets detected. Lower circularity, offset or minimum area.");

    // -------------------------------
    // CONTRAST FILTER (droplet must be brighter than its surround)
    // -------------------------------

    selectImage(originalID);

    for (i = n - 1; i >= 0; i--) {

        roiManager("Select", i);
        inMean = getValue("Mean");
        inArea = getValue("Area");

        run("Enlarge...", "enlarge=4 pixel");
        bigMean = getValue("Mean");
        bigArea = getValue("Area");

        surround = (bigMean * bigArea - inMean * inArea) / (bigArea - inArea);

        if (surround <= 0 || inMean / surround < minContr) {
            roiManager("Select", i);
            roiManager("Delete");
        }
    }

    run("Select None");

    n = roiManager("Count");
    if (n == 0)
        exit("All objects were rejected by the contrast filter. Lower 'Min contrast'.");

    // -------------------------------
    // NAME (AND OPTIONALLY ROUND) THE ROIs
    // -------------------------------

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
    // DISPLAY
    // -------------------------------

    roiManager("Deselect");
    roiManager("Show All without labels");
    roiManager("Set Color", "yellow");
    roiManager("Set Line Width", 1);

    print("Droplets detected: " + n);

    showMessage("Analysis completed",
        "Droplets detected: " + n +
        "\nMean intensity is in the Results table.");
}
