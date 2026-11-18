// Global values returned by ringTest()
var gRingMean = 0;

macro "Automatic Lipid Droplet Circles (circle search v4)" {

    if (nImages == 0)
        exit("Open a fluorescence image first.");

    originalID = getImageID();

    // -------------------------------
    // USER SETTINGS
    // -------------------------------

    Dialog.create("Lipid Droplet Circle Detection v4");

    Dialog.addNumber("Minimum droplet radius (pixels)", 4);
    Dialog.addNumber("Maximum droplet radius (pixels)", 16);
    Dialog.addNumber("Radius step (pixels, 2 = faster)", 1);
    Dialog.addNumber("Ring thickness (pixels)", 3);
    Dialog.addNumber("Smoothing sigma (pixels, 0 = off)", 1);
    Dialog.addNumber("Detection strength (% of strongest ring, lower = more droplets)", 15);
    Dialog.addNumber("Minimum ring completeness (0-1)", 0.75);
    Dialog.addNumber("Minimum spacing (x sum of radii, lower = allow overlap)", 0.70);
    Dialog.addNumber("Minimum ring contrast (% of image bright level, higher = stricter)", 8);
    Dialog.addNumber("Minimum ring brightness (% of image bright level)", 12);
    Dialog.addNumber("Ring must exceed inside/outside by (fraction of ring score, 0-1)", 0.5);
    Dialog.addNumber("ROI radius offset from ring centre-line (pixels, negative = smaller circle)", 0);

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

    if (rMin < 2 || rMax < rMin || rStep < 1 || thickness <= 0 ||
        sigma < 0 || thrPct <= 0 || thrPct >= 100 ||
        minComp < 0 || minComp > 1 || minSpace <= 0 ||
        contrastPct < 0 || brightPct < 0 || margin < 0)
        exit("Check your parameters.");

    ringHalf = thickness / 2;
    gap      = 1;      // space between the ring and its inner/outer comparison bands
    flankW   = 3;      // width of the inner/outer comparison bands
    off      = ringHalf + 2;

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
    // Dim haze / noise is judged against this, so it can never count as a droplet ring.
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
    // response = mean(ring band) - mean(inner + outer comparison bands)
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
        thr = minContrastAbs;      // absolute floor: weak "rings" in haze are ignored

    nc = 0;
    if (maxResp > 0) {
        run("Find Maxima...", "prominence=" + (thr / 2) + " output=[Point Selection]");
        if (selectionType() != -1) {
            getSelectionCoordinates(cx, cy);
            nc = cx.length;
        }
    }

    // Best radius and score for each candidate centre
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
    // - ring must clearly beat BOTH inside and outside in >= minComp of 24 directions
    // - the ring must be genuinely bright (not dim haze)
    // -------------------------------

    selectImage(smID);

    keep = newArray(nc);
    nKeep = 0;

    for (j = 0; j < nc; j++) {

        if (bestScore[j] < thr) continue;

        x0 = cx[j];
        y0 = cy[j];
        rr = bestRad[j];

        // circle must be fully inside the image
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
    // Searches centre +/-2 px and radius +/-2 px (0.5 px steps) for the circle whose
    // average edge brightness is highest = the ring centre-line (not its outer edge).
    // The strict ring test is repeated on the refined circle, so false droplets that
    // only passed because of a slightly misplaced circle are removed.
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

        // second, strict check on the refined circle
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

    setBatchMode(false);
    showProgress(1);
    showStatus("");

    selectImage(originalID);

    if (nFinal == 0)
        exit("No droplets detected. Try lowering 'Minimum ring contrast', " +
             "'Minimum ring brightness' or 'Minimum ring completeness', " +
             "or check the radius range. (If this plane has no droplets, this is correct.)");

    // -------------------------------
    // DRAW CIRCLES AS ROIs
    // Circle goes through the ring centre-line (+ your offset) and is centred on the
    // true pixel centre (+0.5)
    // -------------------------------

    roiManager("Reset");

    for (q = 0; q < nFinal; q++) {
        rr = fr[q] + roiOffset;
        if (rr < 1) rr = 1;
        makeOval(fx[q] + 0.5 - rr, fy[q] + 0.5 - rr, 2 * rr, 2 * rr);
        roiManager("Add");
        roiManager("Select", roiManager("Count") - 1);
        roiManager("Rename", "LD_" + (q + 1));
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

    print("Lipid droplet analysis completed.");
    print("Candidates after first test: " + nAcc);
    print("Detected droplets (after refinement test): " + nFinal);

    showMessage("Analysis completed",
        "Droplets detected: " + nFinal +
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

    innerOut = r - ringHalf - gap;            // inner band: [innerOut - flankW , innerOut)
    outerIn  = r + ringHalf + gap;            // outer band: (outerIn , outerIn + flankW]

    // pass 1: count pixels in each band
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

    // pass 2: write the kernel
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
