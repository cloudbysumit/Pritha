macro "Automatic Lipid Droplet ROI" {

    if (nImages == 0)
        exit("Open a fluorescence image first.");

    original = getTitle();

    // -------------------------------
    // USER SETTINGS
    // -------------------------------

    Dialog.create("Lipid Droplet Detection");

    Dialog.addChoice("Threshold method",
        newArray("Triangle", "Otsu", "Default", "Huang", "Moments"), "Triangle");

    Dialog.addChoice("ROI type",
        newArray("Whole droplet", "Rim only", "Both"), "Whole droplet");

    Dialog.addNumber("Background subtraction radius (pixels, 0 = off)", 25);
    Dialog.addNumber("Gaussian blur sigma (pixels, 0 = off)", 1);
    Dialog.addNumber("Rim ROI thickness (pixels)", 3);
    Dialog.addNumber("Minimum droplet area (pixels)", 15);
    Dialog.addNumber("Maximum droplet area (0 = unlimited)", 0);
    Dialog.addNumber("Minimum circularity (0-1)", 0.40);
    Dialog.addCheckbox("Apply watershed", true);

    Dialog.show();

    method     = Dialog.getChoice();
    roiType    = Dialog.getChoice();
    bgRadius   = Dialog.getNumber();
    blurSigma  = Dialog.getNumber();
    thickness  = Dialog.getNumber();
    minArea    = Dialog.getNumber();
    maxArea    = Dialog.getNumber();
    minCirc    = Dialog.getNumber();
    watershed  = Dialog.getCheckbox();

    if (thickness <= 0 || minArea < 0 ||
        minCirc < 0 || minCirc > 1 ||
        bgRadius < 0 || blurSigma < 0)
        exit("Check your parameters.");

    if (maxArea <= 0)
        maxArea = 999999999;

    half = thickness / 2;

    // -------------------------------
    // SEGMENTATION
    // -------------------------------
    // Duplicates ONLY the currently displayed channel/slice/frame.
    // Make sure the lipid droplet channel is active before running.

    selectWindow(original);
    run("Duplicate...", "title=LD_MASK");
    selectWindow("LD_MASK");

    run("8-bit");

    // Remove uneven cytoplasmic background so droplets stand out
    if (bgRadius > 0)
        run("Subtract Background...", "rolling=" + bgRadius);

    // Light smoothing to reduce noise
    if (blurSigma > 0)
        run("Gaussian Blur...", "sigma=" + blurSigma);

    setOption("BlackBackground", true);
    setAutoThreshold(method + " dark");
    run("Convert to Mask");

    // Fill ring-shaped droplets and close small gaps
    run("Close-");
    run("Fill Holes");

    if (watershed)
        run("Watershed");

    // -------------------------------
    // DETECT DROPLETS (individual ROIs)
    // -------------------------------

    roiManager("Reset");

    run("Analyze Particles...",
        "size=" + minArea + "-" + maxArea + " pixel" +
        " circularity=" + minCirc + "-1.00" +
        " show=Nothing exclude add");

    n = roiManager("Count");

    if (n == 0) {
        selectWindow("LD_MASK");
        close();
        exit("No droplets detected. Adjust threshold, background radius or minimum area.");
    }

    // Name the whole-droplet ROIs
    for (i = 0; i < n; i++) {
        roiManager("Select", i);
        roiManager("Rename", "LD_" + (i + 1));
    }

    // -------------------------------
    // CREATE RIM ROIs (follow the real droplet shape)
    // -------------------------------

    if (roiType != "Whole droplet") {

        for (i = 0; i < n; i++) {

            roiManager("Select", i);
            getSelectionBounds(x, y, w, h);

            // Skip droplets too small for a rim
            if (w <= thickness + 1 || h <= thickness + 1)
                continue;

            // Outer contour
            roiManager("Select", i);
            run("Enlarge...", "enlarge=" + half + " pixel");
            roiManager("Add");
            outerIndex = roiManager("Count") - 1;

            // Inner contour
            roiManager("Select", i);
            run("Enlarge...", "enlarge=-" + half + " pixel");
            roiManager("Add");
            innerIndex = roiManager("Count") - 1;

            // Outer minus inner = ring
            roiManager("Select", newArray(outerIndex, innerIndex));
            roiManager("XOR");
            roiManager("Add");
            ringIndex = roiManager("Count") - 1;
            roiManager("Select", ringIndex);
            roiManager("Rename", "LD_" + (i + 1) + "_Rim");

            // Remove temporary outer and inner ROIs
            roiManager("Select", newArray(outerIndex, innerIndex));
            roiManager("Delete");
        }

        // If only rims wanted, delete the whole-droplet ROIs
        if (roiType == "Rim only") {
            roiManager("Select", Array.getSequence(n));
            roiManager("Delete");
        }
    }

    // -------------------------------
    // REMOVE SEGMENTATION IMAGE
    // -------------------------------

    selectWindow("LD_MASK");
    close();

    // -------------------------------
    // MEASURE ORIGINAL INTENSITY
    // -------------------------------

    selectWindow(original);

    run("Clear Results");
    run("Set Measurements...",
        "area mean min integrated centroid perimeter redirect=None decimal=3");

    finalCount = roiManager("Count");

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
    print("ROI type: " + roiType);
    print("Measured ROIs: " + finalCount);

    showMessage("Analysis completed",
        "Droplets detected: " + n +
        "\nROIs measured: " + finalCount +
        "\nROI type: " + roiType);
}