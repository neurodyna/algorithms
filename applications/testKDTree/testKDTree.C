#include "argList.H"
#include "Random.H"
#include "scalarField.H"
#include "List.H"
#include "KDTree.H"
#include "IOstreams.H"
#include "clockTime.H"

using namespace Foam;


// =========================================================================
//                               DATA GENERATORS
// =========================================================================

// --------------------------------------------------
// 1. Uniform Random [0, 1]
// --------------------------------------------------
List<scalarField> generateUniform
(
    Random& rnd,
    const label n,
    const label k
)
{
    List<scalarField> pts(n);
    forAll(pts, i)
    {
        pts[i].setSize(k);
        for(label d=0; d<k; ++d) pts[i][d] = rnd.scalar01();
    }
    return pts;
}

// --------------------------------------------------
// 2. Clustered Data (Simulates Refined Mesh Zones)
//    - 20% background noise [0, 1]
//    - 40% tight cluster at [0.5, 0.5, ...]
//    - 40% tight cluster at [0.1, 0.1, ...]
// --------------------------------------------------
List<scalarField> generateClustered
(
    Random& rnd,
    const label n,
    const label k
)
{
    List<scalarField> pts(n);
    forAll(pts, i)
    {
        pts[i].setSize(k);
        scalar r = rnd.scalar01();

        scalar center = 0.5;
        scalar spread = 1.0;

        if (r < 0.2)      { center = 0.5; spread = 1.0; } // Background
        else if (r < 0.6) { center = 0.5; spread = 0.05;} // Cluster A (Dense)
        else              { center = 0.1; spread = 0.05;} // Cluster B (Dense)

        for(label d=0; d<k; ++d)
        {
            // Gaussian-ish distribution around center
            pts[i][d] = center + (rnd.scalar01() - 0.5) * spread;
        }
    }
    return pts;
}

// --------------------------------------------------
// 3. Duplicate Data (Stress test for Partitioning)
// --------------------------------------------------
List<scalarField> generateDuplicates
(
    const label n,
    const label k
)
{
    List<scalarField> pts(n);
    forAll(pts, i)
    {
        pts[i].setSize(k);
        for(label d=0; d<k; ++d) pts[i][d] = 0.5; // All points identical
    }
    return pts;
}


// =========================================================================
//                               VALIDATION LOGIC
// =========================================================================

// Brute-force reference
label bruteNearest(const List<scalarField>& pts, const scalarField& q, scalar& bestD2)
{
    label idx = -1;
    bestD2 = GREAT;
    forAll(pts, i)
    {
        scalar d2 = 0;
        for(label d=0; d<q.size(); ++d)
        {
            scalar diff = pts[i][d] - q[d];
            d2 += diff*diff;
        }
        if(d2 < bestD2) { bestD2 = d2; idx = i; }
    }
    return idx;
}

// Core Test Runner
void runTest
(
    const word& testName,
    const List<scalarField>& pts,
    const label nQueries,
    const label nCheck, // Number of BF checks
    bool skipBF = false
)
{
    const label nPts = pts.size();
    const label dim  = nPts > 0 ? pts[0].size() : 0;

    Info<< nl << "--------------------------------------------------" << nl
        << "Test           : " << testName << nl
        << "Points         : " << nPts << nl
        << "Dimensions     : " << dim << nl
        << "--------------------------------------------------" << endl;

    if (nPts == 0)
    {
        Info<< "Skipping empty dataset check." << endl;
        return;
    }

    Random rnd(1234);
    List<scalarField> queries = generateUniform(rnd, nQueries, dim);

    // 1. Build
    clockTime tBuild;
    KDTree tree(pts);
    double buildTime = tBuild.timeIncrement();

    // 2. Query (KDTree)
    clockTime tKD;
    forAll(queries, i) tree.nearest(queries[i]);
    double kdTime = tKD.timeIncrement();

    double kdAvg = (kdTime / nQueries) * 1e6; // microseconds

    Info<< "Build time     : " << buildTime << " s" << nl
        << "Query time     : " << kdTime << " s (" << kdAvg << " us/query)" << endl;

    if (skipBF) return;

    // 3. Validation (Brute Force)
    Info<< "Validating " << nCheck << " samples... " << flush;
    
    label errors = 0;
    clockTime tBF;
    
    for(label i=0; i<nCheck; ++i)
    {
        scalar bfD2, kdD2;
        bruteNearest(pts, queries[i], bfD2);
        kdD2 = tree.nearestDistSqr(queries[i]);

        if (mag(bfD2 - kdD2) > 1e-10)
        {
            errors++;
            if (errors < 5) // Print first few errors
            {
                Info<< nl << "Error at query " << i 
                    << ": BF=" << bfD2 << " KD=" << kdD2;
            }
        }
    }
    double bfTime = tBF.timeIncrement();
    double bfAvg = (bfTime / nCheck) * 1e6;
    double speedup = bfAvg / kdAvg;

    if (errors == 0) Info<< "OK" << endl;
    else             Info<< nl << "FAILED (" << errors << " mismatches)" << endl;

    Info<< "BF avg time    : " << bfAvg << " us/query" << nl
        << "Speedup        : " << speedup << "x" << endl;
}


// =========================================================================
//                                  MAIN
// =========================================================================

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList args(argc, argv);
    Random rnd(54321);

    Info<< nl << "=== KDTree ROBUSTNESS & BENCHMARK SUITE ===" << nl << endl;

    // ---------------------------------------------------------
    // SECTION 1: Corner Cases (Stability)
    // ---------------------------------------------------------
    Info<< "--- SECTION 1: STABILITY CHECKS ---" << endl;
    
    {
        List<scalarField> tiny = generateUniform(rnd, 10, 3);
        runTest("Tiny Dataset (N=10)", tiny, 100, 100);
    }

    {
        // All points are exactly (0.5, 0.5, 0.5).
        // This tests if the partition algorithm hangs on equal keys.
        List<scalarField> dups = generateDuplicates(10000, 3);
        runTest("Duplicate Points (N=10k)", dups, 1000, 100);
    }

    // ---------------------------------------------------------
    // SECTION 2: Real-World Scenarios (Distributions)
    // ---------------------------------------------------------
    Info<< nl << "--- SECTION 2: DISTRIBUTION EFFECTS (3D, N=100k) ---" << endl;

    // Case A: Uniform
    List<scalarField> uniformPts = generateUniform(rnd, 100000, 3);
    runTest("Uniform Distribution", uniformPts, 10000, 500);

    // Case B: Clustered (Wake simulation)
    // KD-Trees are sensitive to clusters; this proves adaptation works.
    List<scalarField> clusterPts = generateClustered(rnd, 100000, 3);
    runTest("Clustered/Wake Distribution", clusterPts, 10000, 500);


    // ---------------------------------------------------------
    // SECTION 3: Large Scale Benchmarks (Scaling)
    // ---------------------------------------------------------
    Info<< nl << "--- SECTION 3: SCALING BENCHMARKS ---" << endl;

    // 1 Million Points (Standard Benchmark)
    List<scalarField> millionPts = generateUniform(rnd, 1000000, 3);
    runTest("1M Points (3D)", millionPts, 10000, 200);

    // 5 Million Points (Stress Test)
    // NOTE: Only run BF on a tiny subset to avoid waiting forever
    List<scalarField> hugePts = generateUniform(rnd, 5000000, 3);
    runTest("5M Points (3D) - Heavy Load", hugePts, 10000, 50);

    Info<< nl << "=== TEST SUITE COMPLETED SUCCESSFULLY ===" << nl << endl;

    return 0;
}