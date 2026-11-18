// Global value returned by ringTest()
var gRingMean = 0;

macro "Automatic Lipid Droplet Circles (circle search v5)" {

    // Open a file chooser if no image is open yet
    if (nImages == 0) {
        path = File.openDialog("Select the lipid droplet image");
        open(path);
    }

    originalID = getImageID();

    // -------------------------------
    // USER SETTINGS
    // -------------------------------

    Dialog.create("Lipid Droplet Circle Detection v5");

    Dialog.addNumber("Minimum droplet radius (pixels)", 4);
    Dialog.addNumber("Maximum droplet radius (pixels)", 16);
    Dialog.addNumber("Radius step (pixels, 2 = faster)", 1);
    Dialog.addNumber("Ring thickness (pixels)", 3);
    Dialog.addNumber("Smoothing sigma (pixels, 0 = off)", 1);
    Dialog.addNumber("Detection strength (% of strongest ring, lower = more droplets)", 15);
    Dialog.addNumber("Minimum ring completeness (0-1)", 0.70);
    Dialog.addNumber("Minimum spacing (x sum of radii, lower = merges more, higher = keeps overlapping droplets)", 0.85);
    Dialog.addNumber("Minimum ring contrast (% of image bright level, higher = stricter)", 4);
    Dialog.addNumber("Minimum ring brightness (% of image bright level)", 4);
    Dialog.addNumber("Ring must exceed inside/outside by (fraction of ring score, 0-1)", 0.2);
    Dialog.addNumber("ROI radius offset from ring centre-line (pixels, negative = smaller circle)", 0);
    Dialog.addChoice("Overlapping droplets",
        newArray("Cut out overlap", "Keep full circles"), "Cut out overlap");
    Dialog.addNumber("Minimum part of a circle that must remain after cut-out (%)", 30);

    Dialog.show();

    rMin        = Dialog.getNumber();
    rMax        = Dialog.getNumber();
    rStep       = Dialog.getNumber();
    thickness   = Dialog.getNumber();
    sigma       = Dialog.getNumber();
    thrPct      = Dialog.getNumber();
    minComp     = Dialog.getNumber();
    minSpace    = Dialog.getNumber();
    contrastPct = Dialog.getNumber();
    brightPct   = Dialog.getNumber();
    margin      = Dialog.getNumber();
    roiOffset   = Dialog.getNumber();
    overlapMode = Dialog.getChoice();
    minRemain   = Dialog.getNumber();

    if (rMin < 2 || rMax < rMin || rStep < 1 || thickness <= 0 ||
        sigma < 0 || thrPct <= 0 || thrPct >= 100 ||
        minComp < 0 || minComp > 1 || minSpace <= 0 ||
        contrastPct < 0 || brightPct < 0 || margin < 0 ||
        minRemain < 0 || minRemain > 100)
        exit("Check your parameters.");

    ringHalf = thickness / 2;
    gap      = 1;      // space between the ring and its inner/outer comparison bands
    flankW   = 3;      // width of the inner/outer comparison bands
    off      = ringHalf + 2;

    setOption("BlackBackground", true);
    setBatchMode(true);

    // -------------------------------
    // PREPARE A SMOOTHED 32-bit COPY OF THE CURRENT PLANE
    // -------------------------------

    selectImage(originalID);

    if (is("hyperstack")) {
        Stack.getPosition(c, z, t);
        run("Duplicate...", "title=LD_SM duplicate channels=" + c +
            " slices=" + z + " frames=" + t);
    } else {
        run("Duplicate...", "title=LD_SM");
    }

    run("Set Scale...", "distance=0 known=0 unit=pixel");   // work in pixels
    run("32-bit");
    if (sigma > 0)
        run("Gaussian Blur...", "sigma=" + sigma);

    smID = getImageID();
    W = getWidth();
    H = getHeight();

    // Reference "bright level" of this plane = 99.5th percentile of the smoothed image.
    getStatistics(area0, mean0, imgMin, imgMax);
    nBins = 1000;
    getHistogram(hv, hc, nBins, imgMin, imgMax);
    cum = 0;
    ref = imgMax;
    for (b = 0; b < nBins; b++) {
        cum = cum + hc[b];
        if (cum >= 0.995 * W * H) {
            ref = hv[b];
            break;
        }
    }
    minContrastAbs = contrastPct / 100 * ref;
    minBrightAbs   = brightPct / 100 * ref;

    // -------------------------------
    // RING RESPONSE FOR EVERY RADIUS
    // -------------------------------

    nR = floor((rMax - rMin) / rStep) + 1;
    respID  = newArray(nR);
    radii   = newArray(nR);
    maxID   = 0;

    for (k = 0; k < nR; k++) {

        r = rMin + k * rStep;
        radii[k] = r;
        showStatus("Searching radius " + r + " px ...");
        showProgress(k, nR);

        kern = buildKernel(r, ringHalf, gap, flankW);

        selectImage(smID);
        run("Duplicate...", "title=LDRESP_" + r);
        respID[k] = getImageID();
        run("Convolve...", "text1=[" + kern + "]");

        if (k == 0) {
            run("Duplicate...", "title=LD_MAX");
            maxID = getImageID();
        } else {
            imageCalculator("Max", maxID, respID[k]);
        }
    }

    // -------------------------------
    // FIND CANDIDATE CENTRES
    // -------------------------------

    selectImage(maxID);
    getStatistics(area, mean, min, maxResp);
    thr = maxResp * thrPct / 100;
    if (thr < minContrastAbs)
        thr = minContrastAbs;

    nc = 0;
    if (maxResp > 0) {
        run("Find Maxima...", "prominence=" + (thr / 2) + " output=[Point Selection]");
        if (selectionType() != -1) {
            getSelectionCoordinates(cx, cy);
            nc = cx.length;
        }
    }

    bestScore = newArray(nc);
    bestRad   = newArray(nc);
    for (j = 0; j < nc; j++)
        bestScore[j] = -1e30;

    for (k = 0; k < nR; k++) {
        selectImage(respID[k]);
        for (j = 0; j < nc; j++) {
            v = getPixel(cx[j], cy[j]);
            if (v > bestScore[j]) {
                bestScore[j] = v;
                bestRad[j]   = radii[k];
            }
        }
    }

    // -------------------------------
    // FIRST RING TEST (on the smoothed copy)
    // -------------------------------

    selectImage(smID);

    keep = newArray(nc);
    nKeep = 0;

    for (j = 0; j < nc; j++) {

        if (bestScore[j] < thr) continue;

        x0 = cx[j];
        y0 = cy[j];
        rr = bestRad[j];

        if (x0 - rr - ringHalf < 0 || y0 - rr - ringHalf < 0 ||
            x0 + rr + ringHalf >= W || y0 + rr + ringHalf >= H)
            continue;

        frac = ringTest(x0, y0, rr, off, margin * bestScore[j]);

        if (frac >= minComp && gRingMean >= minBrightAbs) {
            keep[nKeep] = j;
            nKeep++;
        }
    }

    // -------------------------------
    // REMOVE DUPLICATES (strongest first)
    // Only near-identical detections of the SAME droplet are merged.
    // Partly overlapping neighbouring droplets are kept.
    // -------------------------------

    keepScore = newArray(nKeep);
    for (i = 0; i < nKeep; i++)
        keepScore[i] = bestScore[keep[i]];

    order = newArray(0);
    if (nKeep > 0)
        order = Array.rankPositions(keepScore);   // ascending
    accepted = newArray(nKeep);
    nAcc = 0;

    for (o = nKeep - 1; o >= 0; o--) {

        j = keep[order[o]];
        ok = true;

        for (q = 0; q < nAcc; q++) {
            jj = accepted[q];
            d = sqrt((cx[j] - cx[jj]) * (cx[j] - cx[jj]) +
                     (cy[j] - cy[jj]) * (cy[j] - cy[jj]));
            if (d < minSpace * (bestRad[j] + bestRad[jj])) {
                ok = false;
                break;
            }
        }

        if (ok) {
            accepted[nAcc] = j;
            nAcc++;
        }
    }

    // -------------------------------
    // REFINE EACH CIRCLE ONTO THE BRIGHT RING, THEN TEST AGAIN
    // -------------------------------

    selectImage(smID);

    fx = newArray(nAcc);
    fy = newArray(nAcc);
    fr = newArray(nAcc);
    nFinal = 0;
    nAng = 36;
    PI2 = 2 * PI;

    for (q = 0; q < nAcc; q++) {

        showStatus("Refining circle " + (q + 1) + " of " + nAcc);

        j  = accepted[q];
        bx = cx[j];
        by = cy[j];
        br = bestRad[j];
        bestProf = -1e30;
        px0 = bx;
        py0 = by;
        rc0 = br;

        for (dy = -2; dy <= 2; dy++) {
            for (dx = -2; dx <= 2; dx++) {
                for (rc = br - 2; rc <= br + 2; rc = rc + 0.5) {

                    if (rc < 2) continue;
                    px = bx + dx;
                    py = by + dy;

                    if (px - rc - ringHalf < 0 || py - rc - ringHalf < 0 ||
                        px + rc + ringHalf >= W || py + rc + ringHalf >= H)
                        continue;

                    ps = 0;
                    for (a = 0; a < nAng; a++) {
                        ang = a * PI2 / nAng;
                        ps = ps + getPixel(px + rc * cos(ang), py + rc * sin(ang));
                    }

                    if (ps > bestProf) {
                        bestProf = ps;
                        px0 = px;
                        py0 = py;
                        rc0 = rc;
                    }
                }
            }
        }

        frac = ringTest(px0, py0, rc0, off, margin * bestScore[j]);

        if (frac >= minComp && gRingMean >= minBrightAbs) {
            fx[nFinal] = px0;
            fy[nFinal] = py0;
            fr[nFinal] = rc0;
            nFinal++;
        }
    }

    // -------------------------------
    // CLOSE TEMPORARY IMAGES
    // -------------------------------

    for (k = 0; k < nR; k++) {
        selectImage(respID[k]);
        close();
    }
    selectImage(maxID);
    close();
    selectImage(smID);
    close();

    if (nFinal == 0) {
        setBatchMode(false);
        showStatus("");
        selectImage(originalID);
        exit("No droplets detected. Try lowering 'Minimum ring contrast', " +
             "'Minimum ring brightness' or 'Minimum ring completeness', " +
             "or check the radius range. (If this plane has no droplets, this is correct.)");
    }

    // -------------------------------
    // BUILD THE ROIs
    // "Keep full circles": one circle per droplet.
    // "Cut out overlap":   the area shared by two or more circles is removed from each,
    //                      so every droplet keeps only its own area (set subtraction).
    // -------------------------------

    roiManager("Reset");
    nOverlap = 0;
    nDropped = 0;

    if (overlapMode == "Keep full circles") {

        for (q = 0; q < nFinal; q++) {
            rr = fr[q] + roiOffset;
            if (rr < 1) rr = 1;
            selectImage(originalID);
            makeOval(fx[q] + 0.5 - rr, fy[q] + 0.5 - rr, 2 * rr, 2 * rr);
            roiManager("Add");
            roiManager("Select", roiManager("Count") - 1);
            roiManager("Rename", "LD_" + (q + 1));
        }

    } else {

        // 1) count how many circles cover each pixel
        newImage("LD_CNT", "8-bit black", W, H, 1);
        cntID = getImageID();

        for (q = 0; q < nFinal; q++) {
            rr = fr[q] + roiOffset;
            if (rr < 1) rr = 1;
            makeOval(fx[q] + 0.5 - rr, fy[q] + 0.5 - rr, 2 * rr, 2 * rr);
            run("Add...", "value=1");
        }
        run("Select None");

        // 2) overlap mask = pixels covered by 2 or more circles (255), everything else 0
        run("Subtract...", "value=1");
        run("Multiply...", "value=255");
        rename("LD_OV");

        // 3) for each droplet: its full circle minus the overlap mask
        newImage("LD_TMP", "8-bit black", W, H, 1);
        tmpID = getImageID();

        for (q = 0; q < nFinal; q++) {

            rr = fr[q] + roiOffset;
            if (rr < 1) rr = 1;

            selectImage(tmpID);
            run("Select All");
            setColor(0);
            fill();
            makeOval(fx[q] + 0.5 - rr, fy[q] + 0.5 - rr, 2 * rr, 2 * rr);
            setColor(255);
            fill();
            run("Select None");

            imageCalculator("Subtract", "LD_TMP", "LD_OV");

            selectImage(tmpID);
            setThreshold(255, 255);
            run("Create Selection");

            if (selectionType() == -1) {
                nDropped++;
                continue;
            }

            remain = getValue("Area");
            full   = PI * rr * rr;

            if (remain < minRemain / 100 * full) {
                nDropped++;
                run("Select None");
                continue;
            }

            if (remain < 0.98 * full)
                nOverlap++;

            roiManager("Add");
            roiManager("Select", roiManager("Count") - 1);
            roiManager("Rename", "LD_" + (q + 1));
            run("Select None");
        }

        selectImage(tmpID);
        close();
        selectImage(cntID);
        close();
    }

    setBatchMode(false);
    showProgress(1);
    showStatus("");

    finalCount = roiManager("Count");

    if (finalCount == 0)
        exit("All droplets were removed by the overlap cut-out. Lower 'Minimum part of a circle that must remain' " +
             "or choose 'Keep full circles'.");

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

    print("Lipid droplet analysis completed.");
    print("Droplets detected: " + nFinal);
    print("ROIs measured: " + finalCount);
    print("ROIs with overlap cut out: " + nOverlap);
    print("Droplets dropped (too little area left): " + nDropped);

    showMessage("Analysis completed",
        "Droplets detected: " + nFinal +
        "\nROIs measured: " + finalCount +
        "\nROIs with overlap cut out: " + nOverlap +
        "\nMean intensity is in the Results table.");
}

// ------------------------------------------------------------
// Ring test on the currently selected (smoothed) image.
// Returns the fraction of 24 directions where the ring beats BOTH the inside and
// the outside by more than 'minGap'. Sets gRingMean = average ring brightness.
// ------------------------------------------------------------
function ringTest(x0, y0, rr, off, minGap) {

    PI2 = 2 * PI;
    good = 0;
    ringSum = 0;

    for (a = 0; a < 24; a++) {
        ang = a * PI2 / 24;
        ca = cos(ang);
        sa = sin(ang);

        ringV = getPixel(x0 + rr * ca, y0 + rr * sa);
        v2 = getPixel(x0 + (rr - 1) * ca, y0 + (rr - 1) * sa);
        v3 = getPixel(x0 + (rr + 1) * ca, y0 + (rr + 1) * sa);
        if (v2 > ringV) ringV = v2;
        if (v3 > ringV) ringV = v3;

        inV  = getPixel(x0 + (rr - off) * ca, y0 + (rr - off) * sa);
        outV = getPixel(x0 + (rr + off) * ca, y0 + (rr + off) * sa);

        ringSum = ringSum + ringV;
        sideMax = inV;
        if (outV > sideMax) sideMax = outV;

        if (ringV > sideMax + minGap)
            good++;
    }

    gRingMean = ringSum / 24;
    return good / 24;
}

// ------------------------------------------------------------
// Builds a zero-sum ring-detector kernel as text for Process > Filters > Convolve
//   +1/N on the ring band (radius r +/- ringHalf)
//   -1/M on the inner and outer comparison bands
// ------------------------------------------------------------
function buildKernel(r, ringHalf, gap, flankW) {

    half = Math.ceil(r + ringHalf + gap + flankW);

    innerOut = r - ringHalf - gap;
    outerIn  = r + ringHalf + gap;

    nRing = 0;
    nFlank = 0;
    for (y = -half; y <= half; y++) {
        for (x = -half; x <= half; x++) {
            d = sqrt(x * x + y * y);
            if (abs(d - r) <= ringHalf)
                nRing++;
            else if ((d >= innerOut - flankW && d < innerOut) ||
                     (d > outerIn && d <= outerIn + flankW))
                nFlank++;
        }
    }

    wRing  = 1 / nRing;
    wFlank = -1 / nFlank;

    s = "";
    for (y = -half; y <= half; y++) {
        row = "";
        for (x = -half; x <= half; x++) {
            d = sqrt(x * x + y * y);
            w = 0;
            if (abs(d - r) <= ringHalf)
                w = wRing;
            else if ((d >= innerOut - flankW && d < innerOut) ||
                     (d > outerIn && d <= outerIn + flankW))
                w = wFlank;
            row = row + d2s(w, 6) + " ";
        }
        s = s + row + "\n";
    }
    return s;
}
