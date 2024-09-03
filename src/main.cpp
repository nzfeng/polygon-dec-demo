#include "geometrycentral/surface/exact_geodesics.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/utilities.h"

#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "signed_heat_polygon.h"

#include "args/args.hxx"
#include "imgui.h"

#include <chrono>
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// == Polyscope data
polyscope::SurfaceMesh* psMesh;
std::vector<Vector3> VBASISX, VBASISY;

// == other data
double MAX_DIAGONAL_LENGTH = 0.;
double MEAN_EDGE_LENGTH = 0.;
std::vector<std::vector<Vertex>> CURVES;
std::unique_ptr<SignedHeatPolygon> signedHeatSolver;

std::vector<std::vector<Vertex>> readCurves(const std::string& filename) {

    std::vector<std::vector<Vertex>> curves;
    std::ifstream curr_file(filename.c_str());
    std::string line;
    std::string X;
    size_t idx;
    bool flip = false;
    bool read = false;
    std::vector<bool> curveFlips;

    if (curr_file.is_open()) {
        while (!curr_file.eof()) {
            getline(curr_file, line);
            // Ignore any newlines
            if (line == "") {
                continue;
            }
            std::istringstream iss(line);
            iss >> X;
            if (X == "signed_curve") {
                iss >> flip;
                curveFlips.push_back(flip);
                curves.emplace_back();
                read = true;
            } else if (X == "unsigned_curve") {
                read = false;
            } else if (X == "unsigned_point") {
                read = false;
            }
            if (read) {
                if (X == "v") {
                    iss >> idx;
                    curves.back().emplace_back(mesh->vertex(idx));
                } else if (X == "l") {
                    while (true) {
                        if (iss.eof()) break;
                        iss >> idx;
                        idx -= 1; // OBJ elements are 1-indexed, whereas geometry-central uses 0-indexing
                        curves.back().emplace_back(mesh->vertex(idx));
                    }
                }
            }
        }
        curr_file.close();

        for (size_t i = 0; i < curves.size(); i++) {
            if (curveFlips[i]) {
                std::reverse(curves[i].begin(), curves[i].end());
            }
        }

    } else {
        std::cerr << "Could not open file <" << filename << ">." << std::endl;
    }
    return curves;
}

void displayInput(const std::vector<std::vector<Vertex>>& curves, const std::string& name = "input",
                  bool display = true) {

    std::vector<Vector3> nodes;
    std::vector<std::array<size_t, 2>> edges;

    size_t offset = 0;
    for (const auto& curve : curves) {
        size_t N = curve.size();
        for (size_t i = 0; i < N - 1; i++) {
            nodes.push_back(geometry->vertexPositions[curve[i]]);
            edges.push_back({offset + i, offset + i + 1});
        }
        nodes.push_back(geometry->vertexPositions[curve[N - 1]]);
        offset += N;
    }

    auto psCurve = polyscope::registerCurveNetwork(name, nodes, edges);
    psCurve->setEnabled(display);
    psCurve->setColor({1, 0.3, 0});
}

void timing() {

    std::cerr << "timing()..." << std::endl;
    std::chrono::time_point<high_resolution_clock> t1, t2;
    std::chrono::milliseconds ms_int;

    std::cerr << "\n";
    t1 = high_resolution_clock::now();
    geometry->requireCotanLaplacian();
    t2 = high_resolution_clock::now();
    ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cerr << "\trequireCotanLaplacian(): " << ms_int.count() << "ms" << std::endl;

    t1 = high_resolution_clock::now();
    geometry->requireVertexLumpedMassMatrix();
    t2 = high_resolution_clock::now();
    ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cerr << "\trequireVertexLumpedMassMatrix(): " << ms_int.count() << "ms" << std::endl;

    t1 = high_resolution_clock::now();
    geometry->requireDECOperators();
    t2 = high_resolution_clock::now();
    ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cerr << "\trequireDECOperators(): " << ms_int.count() << "ms" << std::endl;

    std::cerr << "\n";

    t1 = high_resolution_clock::now();
    geometry->requireSimplePolygonLaplacian();
    t2 = high_resolution_clock::now();
    ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cerr << "\trequireSimplePolygonLaplacian(): " << ms_int.count() << "ms" << std::endl;

    t1 = high_resolution_clock::now();
    geometry->requireSimplePolygonVertexLumpedMassMatrix();
    t2 = high_resolution_clock::now();
    ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cerr << "\trequireSimplePolygonVertexLumpedMassMatrix(): " << ms_int.count() << "ms" << std::endl;

    t1 = high_resolution_clock::now();
    geometry->requireSimplePolygonVertexGalerkinMassMatrix();
    t2 = high_resolution_clock::now();
    ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cerr << "\trequireSimplePolygonVertexGalerkinMassMatrix(): " << ms_int.count() << "ms" << std::endl;

    std::cerr << "\n";

    t1 = high_resolution_clock::now();
    geometry->requirePolygonLaplacian();
    t2 = high_resolution_clock::now();
    ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cerr << "\trequirePolygonLaplacian(): " << ms_int.count() << "ms" << std::endl;

    t1 = high_resolution_clock::now();
    geometry->requirePolygonGradientMatrix();
    t2 = high_resolution_clock::now();
    ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cerr << "\trequirePolygonGradientMatrix(): " << ms_int.count() << "ms" << std::endl;

    t1 = high_resolution_clock::now();
    geometry->requirePolygonDivergenceMatrix();
    t2 = high_resolution_clock::now();
    ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cerr << "\trequirePolygonDivergenceMatrix(): " << ms_int.count() << "ms" << std::endl;

    t1 = high_resolution_clock::now();
    geometry->requirePolygonVertexLumpedMassMatrix();
    t2 = high_resolution_clock::now();
    ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cerr << "\trequirePolygonVertexLumpedMassMatrix(): " << ms_int.count() << "ms" << std::endl;

    t1 = high_resolution_clock::now();
    geometry->requirePolygonVertexConnectionLaplacian();
    t2 = high_resolution_clock::now();
    ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cerr << "\trequirePolygonVertexConnectionLaplacian(): " << ms_int.count() << "ms" << std::endl;

    t1 = high_resolution_clock::now();
    geometry->requirePolygonDECOperators();
    t2 = high_resolution_clock::now();
    ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cerr << "\trequirePolygonDECOperators(): " << ms_int.count() << "ms" << std::endl;

    std::cerr << "\tDone testing." << std::endl;
}

/* Test that the lumped mass matrices correspond with their unlumped versions. */
void testMassLumping() {
    std::cerr << "testMassLumping()..." << std::endl;

    geometry->requireSimplePolygonVertexLumpedMassMatrix();
    geometry->requireSimplePolygonVertexGalerkinMassMatrix();

    double epsilon = 1e-8;
    SparseMatrix<double> M_T = geometry->simplePolygonVertexGalerkinMassMatrix.transpose();
    SparseMatrix<double> M = geometry->simplePolygonVertexLumpedMassMatrix;
    for (size_t i = 0; i < mesh->nVertices(); i++) {
        double rowSum = 0.;
        for (SparseMatrix<double>::InnerIterator it(M_T, i); it; ++it) {
            rowSum += it.value();
        }
        assert(std::abs(rowSum - M.coeffRef(i, i)) < epsilon);
    }

    geometry->unrequireSimplePolygonVertexLumpedMassMatrix();
    geometry->unrequireSimplePolygonVertexGalerkinMassMatrix();
    std::cerr << "\tDone testing." << std::endl;
}

/* Check that polygon Laplacians and mass matrices coincide with simplicial versions on a triangle mesh. */
void testOnTriangleMeshes() {

    geometry->requireVertexGalerkinMassMatrix();
    geometry->requireVertexLumpedMassMatrix();
    geometry->requireCotanLaplacian();
    geometry->requireSimplePolygonLaplacian();
    geometry->requireSimplePolygonVertexLumpedMassMatrix();
    geometry->requireSimplePolygonVertexGalerkinMassMatrix();
    geometry->requirePolygonLaplacian();
    geometry->requirePolygonVertexLumpedMassMatrix();

    SparseMatrix<double>& L = geometry->cotanLaplacian;
    SparseMatrix<double>& L_VR = geometry->simplePolygonLaplacian;
    SparseMatrix<double>& L_VEM = geometry->polygonLaplacian;
    SparseMatrix<double>& M = geometry->vertexGalerkinMassMatrix;
    SparseMatrix<double>& ML = geometry->vertexLumpedMassMatrix;
    SparseMatrix<double>& M_VR = geometry->simplePolygonVertexGalerkinMassMatrix;
    SparseMatrix<double>& M_VRL = geometry->simplePolygonVertexLumpedMassMatrix;
    SparseMatrix<double>& M_VEM = geometry->polygonVertexLumpedMassMatrix;

    std::cerr << "Laplace matrices:" << std::endl;
    std::cerr << "\t|L - L_VR|: " << (L - L_VR).norm() << "\t" << (L - L_VR).norm() / L.norm() << std::endl;
    std::cerr << "\t|L - L_VEM|: " << (L - L_VEM).norm() << "\t" << (L - L_VEM).norm() / L.norm() << std::endl;
    std::cerr << "Mass matrices:" << std::endl;
    std::cerr << "\t|M - M_VR|: " << (M - M_VR).norm() << "\t" << (M - M_VR).norm() / M.norm() << std::endl;
    std::cerr << "\t|ML - M_VRL|: " << (ML - M_VRL).norm() << "\t" << (ML - M_VRL).norm() / ML.norm() << std::endl;
    std::cerr << "\t|ML - M_VEM|: " << (ML - M_VEM).norm() << "\t" << (ML - M_VEM).norm() / ML.norm() << std::endl;
    std::cerr << "Unlumped vs. lumped matrices (should have significant but identical errors):" << std::endl;
    std::cerr << "\t|M - M_VRL|: " << (M - M_VRL).norm() << "\t" << (M - M_VRL).norm() / M.norm() << std::endl;
    std::cerr << "\t|M - M_VEM|: " << (M - M_VEM).norm() << "\t" << (M - M_VEM).norm() / M.norm() << std::endl;
    std::cerr << "\t|M - ML|: " << (M - ML).norm() << "\t" << (M - ML).norm() / M.norm() << std::endl;

    geometry->unrequireVertexGalerkinMassMatrix();
    geometry->unrequireVertexLumpedMassMatrix();
    geometry->unrequireCotanLaplacian();
    geometry->unrequireSimplePolygonLaplacian();
    geometry->unrequireSimplePolygonVertexLumpedMassMatrix();
    geometry->unrequireSimplePolygonVertexGalerkinMassMatrix();
    geometry->unrequirePolygonLaplacian();
    geometry->unrequirePolygonVertexLumpedMassMatrix();
}

/* Check that gradient, divergence, Laplacian can be assembled as expected from DEC operators. */
void testDECOperators() {
    std::cerr << "testDECOperators()..." << std::endl;
    geometry->requirePolygonDECOperators();
    geometry->requirePolygonLaplacian();

    double epsilon = 1e-8;
    SparseMatrix<double>& L = geometry->polygonLaplacian;
    SparseMatrix<double>& h0 = geometry->polygonHodge0;
    SparseMatrix<double>& h0Inv = geometry->polygonHodge0Inverse;
    SparseMatrix<double>& h1 = geometry->polygonHodge1;
    SparseMatrix<double>& h2 = geometry->polygonHodge2;
    SparseMatrix<double>& h2Inv = geometry->polygonHodge2Inverse;
    SparseMatrix<double>& d0 = geometry->polygonD0;
    SparseMatrix<double>& d1 = geometry->polygonD1;
    assert((L - d0.transpose() * h1 * d0).norm() < epsilon);

    if (mesh->isTriangular()) {
        geometry->requireDECOperators();
        geometry->requireCotanLaplacian();

        std::cerr << "\t|C - L|: " << (geometry->cotanLaplacian - L).norm() << std::endl;
        std::cerr << "\t|C - d*d|: "
                  << (geometry->cotanLaplacian - geometry->d0.transpose() * geometry->hodge1 * geometry->d0).norm()
                  << std::endl;
        std::cerr << "\t|h0 - ph0|: " << (h0 - geometry->hodge0).norm() << std::endl;
        std::cerr << "\t|h2 - ph2|: " << (h2 - geometry->hodge2).norm() << std::endl;

        geometry->unrequireDECOperators();
        geometry->unrequireCotanLaplacian();
    }

    // TODO: Check gradient & divergence with sharp and flat operators?

    geometry->unrequirePolygonDECOperators();
    geometry->unrequirePolygonLaplacian();
    std::cerr << "\tDone testing." << std::endl;
}

/* Solve Poisson problems using both sets of polygon operators. */
void solvePoissonProblems() {
    std::cerr << "solvePoissonProblems()..." << std::endl;
    // Randomly choose some vertex sources.
    size_t V = mesh->nVertices();
    int nSources = randomInt(1, 5);
    std::vector<Vector3> sourcePositions(nSources);
    Vector<double> rho = Vector<double>::Zero(V);
    for (int i = 0; i < nSources; i++) {
        size_t vIdx = randomIndex(V);
        sourcePositions[i] = geometry->vertexPositions[vIdx];
        double weight = randomInt(1, 5);
        weight *= (unitRand() < 0.5) ? 1. : -1.;
        rho[vIdx] += weight;
    }

    // Plot sources.
    polyscope::registerPointCloud("Poisson sources", sourcePositions);

    geometry->requireSimplePolygonLaplacian();
    geometry->requireSimplePolygonVertexLumpedMassMatrix();
    geometry->requireSimplePolygonVertexGalerkinMassMatrix();
    geometry->requirePolygonLaplacian();
    geometry->requirePolygonVertexLumpedMassMatrix();

    SparseMatrix<double> L_VR = geometry->simplePolygonLaplacian;
    SparseMatrix<double> L_VEM = geometry->polygonLaplacian;
    SparseMatrix<double> M_VRL = geometry->simplePolygonVertexLumpedMassMatrix;
    SparseMatrix<double> M_VR = geometry->simplePolygonVertexGalerkinMassMatrix;
    SparseMatrix<double> M_VEM = geometry->polygonVertexLumpedMassMatrix;
    double totalRho_VR = (M_VR * rho).sum();
    double totalRho_VRL = (M_VRL * rho).sum();
    double totalRho_VEM = (M_VEM * rho).sum();
    Vector<double> ones = Vector<double>::Ones(V);
    double totalArea_VR = (M_VR * ones).sum();
    double totalArea_VEM = (M_VEM * ones).sum();
    std::cerr << "\ttotalArea_VR: " << totalArea_VR << std::endl;
    std::cerr << "\ttotalArea_VEM: " << totalArea_VEM << std::endl;
    Vector<double> rhoBar_VR = Vector<double>::Ones(V) * (totalRho_VR / totalArea_VR);
    Vector<double> rhs_VR = M_VR * (rhoBar_VR - rho);
    Vector<double> rhoBar_VRL = Vector<double>::Ones(V) * (totalRho_VRL / totalArea_VR);
    Vector<double> rhs_VRL = M_VRL * (rhoBar_VRL - rho);
    Vector<double> rhoBar_VEM = Vector<double>::Ones(V) * (totalRho_VEM / totalArea_VEM);
    Vector<double> rhs_VEM = M_VEM * (rhoBar_VEM - rho);
    shiftDiagonal(L_VR, 1e-8);
    shiftDiagonal(L_VEM, 1e-8);
    Vector<double> VRSoln = solvePositiveDefinite(L_VR, rhs_VR);
    Vector<double> VRLSoln = solvePositiveDefinite(L_VR, rhs_VRL);
    Vector<double> VEMSoln = solvePositiveDefinite(L_VEM, rhs_VEM);
    std::cerr << "\t[poisson] |VRSoln - VRLSoln| / |VRSoln|: " << (VRSoln - VRLSoln).norm() / VRSoln.norm()
              << std::endl;
    std::cerr << "\t[poisson] |VRSoln - VRLSoln| / |VRLSoln|: " << (VRSoln - VRLSoln).norm() / VRLSoln.norm()
              << std::endl;
    std::cerr << "\t[poisson] |VRSoln - VEMSoln| / |VRSoln|: " << (VRSoln - VEMSoln).norm() / VRSoln.norm()
              << std::endl;
    std::cerr << "\t[poisson] |VRSoln - VEMSoln| / |VEMSoln|: " << (VRSoln - VEMSoln).norm() / VEMSoln.norm()
              << std::endl;
    std::cerr << "\n";

    if (mesh->isTriangular()) {
        // Solve using standard cotan operator.
        geometry->requireVertexLumpedMassMatrix();
        geometry->requireCotanLaplacian();
        geometry->requireFaceAreas();

        SparseMatrix<double> C = geometry->cotanLaplacian;
        SparseMatrix<double> M = geometry->vertexLumpedMassMatrix;
        double totalArea = 0.;
        for (Face f : mesh->faces()) totalArea += geometry->faceAreas[f];
        std::cerr << "\ttotalArea: " << totalArea << std::endl;
        double totalRho = (M * rho).sum();
        Vector<double> rhoBar = Vector<double>::Ones(V) * (totalRho / totalArea);
        Vector<double> rhs = M * (rhoBar - rho);
        shiftDiagonal(C, 1e-8);
        Vector<double> triSoln = solvePositiveDefinite(C, rhs);
        // On triangle meshes, both polygon Laplacians and mass matrices (and hence solutions) should coincide with the
        // cotan solution.
        std::cerr << "\ttotalRho: " << totalRho << "\tVR: " << totalRho_VR << "\tVRL: " << totalRho_VRL
                  << "\tVEM: " << totalRho_VEM << std::endl;
        std::cerr << "\t|rhs - rhs_VR|: " << (rhs - rhs_VR).norm() << "\t" << (rhs - rhs_VR).norm() / rhs.norm()
                  << std::endl;
        std::cerr << "\t|rhs - rhs_VRL|: " << (rhs - rhs_VRL).norm() << "\t" << (rhs - rhs_VRL).norm() / rhs.norm()
                  << std::endl;
        std::cerr << "\t|rhs - rhs_VEM|: " << (rhs - rhs_VEM).norm() << "\t" << (rhs - rhs_VEM).norm() / rhs.norm()
                  << std::endl;
        psMesh->addVertexScalarQuantity("tri", triSoln);
        psMesh->addVertexScalarQuantity("VR", VRSoln);
        psMesh->addVertexScalarQuantity("VR (lumped)", VRLSoln);
        psMesh->addVertexScalarQuantity("VEM", VEMSoln);
        std::cerr << "\t[poisson] |VR - truth|: " << (triSoln - VRSoln).norm() << "\t"
                  << (triSoln - VRSoln).norm() / triSoln.norm() << std::endl;
        std::cerr << "\t[poisson] |VR (lumped) - truth|: " << (triSoln - VRLSoln).norm() << "\t"
                  << (triSoln - VRLSoln).norm() / triSoln.norm() << std::endl;
        std::cerr << "\t[poisson] |VEM - truth|: " << (triSoln - VEMSoln).norm() << "\t"
                  << (triSoln - VEMSoln).norm() / triSoln.norm() << std::endl;

        geometry->unrequireCotanLaplacian();
        geometry->unrequireVertexLumpedMassMatrix();
        geometry->unrequireFaceAreas();
    }

    geometry->unrequireSimplePolygonLaplacian();
    geometry->unrequireSimplePolygonVertexLumpedMassMatrix();
    geometry->unrequirePolygonLaplacian();
    geometry->unrequirePolygonVertexLumpedMassMatrix();
    std::cerr << "\tDone testing." << std::endl;
}

/* Implement the scalar heat method -- good way to test the gradient, divergence, and Laplacian. */
void solveGeodesicDistance() {
    std::cerr << "solveGeodesicDistance()..." << std::endl;
    // Randomly choose some vertex sources.
    size_t V = mesh->nVertices();
    size_t F = mesh->nFaces();
    Vector<double> rho = Vector<double>::Zero(V);
    const int nSources = randomInt(1, 5);
    std::vector<Vector3> sourcePositions(nSources);
    std::vector<Vertex> sources(nSources);
    for (int i = 0; i < nSources; i++) {
        size_t idx = randomIndex(V);
        sourcePositions[i] = geometry->vertexPositions[idx];
        sources[i] = mesh->vertex(idx);
        rho[idx] += 1.;
    }

    // Plot sources.
    polyscope::registerPointCloud("geodesic sources", sourcePositions);

    geometry->requireVertexIndices();
    geometry->requirePolygonVertexLumpedMassMatrix();
    geometry->requirePolygonLaplacian();
    geometry->requirePolygonGradientMatrix();
    geometry->requirePolygonDivergenceMatrix();
    geometry->requirePolygonDECOperators();

    // Flow heat.
    SparseMatrix<double> M = geometry->polygonVertexLumpedMassMatrix;
    SparseMatrix<double> L = geometry->polygonLaplacian;
    double shortTime = MAX_DIAGONAL_LENGTH * MAX_DIAGONAL_LENGTH;
    SparseMatrix<double> LHS = M + shortTime * L;
    Vector<double> X = solvePositiveDefinite(LHS, rho);

    // Normalize gradient.
    Vector<double> Y = geometry->polygonGradientMatrix * X; // 3|F|
    for (size_t i = 0; i < F; i++) {
        Vector3 g = {Y[3 * i], Y[3 * i + 1], Y[3 * i + 2]};
        g /= g.norm();
        Y[3 * i] = g[0];
        Y[3 * i + 1] = g[1];
        Y[3 * i + 2] = g[2];
    }

    // Integrate.
    Vector<double> div = geometry->polygonDivergenceMatrix * Y;
    shiftDiagonal(L, 1e-8);
    Vector<double> VEMDistances = solvePositiveDefinite(L, div);
    VEMDistances *= -1.; // since div * grad gives positive Laplacian

    // Shift solution.
    double avgVEM = 0.;
    for (const Vertex& v : sources) {
        size_t vIdx = geometry->vertexIndices[v];
        avgVEM += VEMDistances[vIdx];
    }
    avgVEM /= nSources;
    VEMDistances -= avgVEM * Vector<double>::Ones(V);

    geometry->unrequireVertexIndices();
    geometry->unrequirePolygonLaplacian();
    geometry->unrequirePolygonGradientMatrix();
    geometry->unrequirePolygonDivergenceMatrix();
    geometry->unrequirePolygonVertexLumpedMassMatrix();
    geometry->unrequirePolygonDECOperators();

    psMesh->addVertexSignedDistanceQuantity("VEM", VEMDistances);

    if (mesh->isTriangular()) {
        // Solve using standard operators.
        HeatMethodDistanceSolver triSolver(*geometry);
        VertexData<double> hmDistances = triSolver.computeDistance(sources);

        GeodesicAlgorithmExact mmp(*mesh, *geometry);
        mmp.propagate(sources);
        VertexData<double> mmpDistances = mmp.getDistanceFunction();

        psMesh->addVertexSignedDistanceQuantity("heat method", hmDistances);
        std::cerr << "\t[distance] |VEM - truth|: " << (mmpDistances.toVector() - VEMDistances).norm() << "\t"
                  << (mmpDistances.toVector() - VEMDistances).norm() / mmpDistances.toVector().norm() << std::endl;
        std::cerr << "\t[distance] |heat method - truth|: " << (mmpDistances.toVector() - hmDistances.toVector()).norm()
                  << "\t" << (mmpDistances.toVector() - hmDistances.toVector()).norm() / mmpDistances.toVector().norm()
                  << std::endl;
    }
    std::cerr << "\tDone testing." << std::endl;
}

void solveVectorHeatMethod() {

    // Randomly choose some vertex sources with random magnitudes.
    size_t V = mesh->nVertices();
    const int nSources = randomInt(2, 5);
    std::vector<std::tuple<Vertex, double>> sourceMagnitudes(nSources);
    std::vector<std::tuple<Vertex, Vector2>> sourceVectors(nSources);
    std::vector<Vector3> sourcePositions(nSources);
    VertexData<Vector3> vectorSources(*mesh);
    for (int i = 0; i < nSources; i++) {
        size_t idx = randomIndex(V);
        Vertex v = mesh->vertex(idx);
        sourceMagnitudes[i] = std::make_tuple(v, randomReal(1., 5.));
        Vector2 vec = {randomReal(-1., 1.), randomReal(-1., 1.)};
        vec /= vec.norm();
        sourceVectors[i] = std::make_tuple(v, vec);
        sourcePositions[i] = geometry->vertexPositions[idx];
        vectorSources[v] = vec[0] * VBASISX[idx] + vec[1] * VBASISX[idx];
    }

    // Plot sources.
    polyscope::registerPointCloud("source locations", sourcePositions);

    geometry->requireVertexIndices();
    geometry->requirePolygonLaplacian();
    geometry->requirePolygonVertexConnectionLaplacian();
    geometry->requirePolygonVertexLumpedMassMatrix();

    // Encode initial data
    Vector<double> rho_magnitudes = Vector<double>::Zero(V);
    Vector<double> rho = Vector<double>::Zero(V);
    Vector<std::complex<double>> rho_vectors = Vector<std::complex<double>>::Zero(V);
    for (int i = 0; i < nSources; i++) {
        size_t vIdx = geometry->vertexIndices[std::get<0>(sourceMagnitudes[i])];
        rho[vIdx] += 1.0;
        rho_magnitudes[vIdx] += std::get<1>(sourceMagnitudes[i]);
        rho_vectors[vIdx] += std::complex<double>(std::get<1>(sourceVectors[i]));
    }

    // Scalar extension
    double shortTime = MAX_DIAGONAL_LENGTH * MAX_DIAGONAL_LENGTH;
    SparseMatrix<double> M = geometry->polygonVertexLumpedMassMatrix;
    SparseMatrix<double> L = geometry->polygonLaplacian;
    SparseMatrix<double> heatOp = M + shortTime * L;
    PositiveDefiniteSolver<double> heatSolver(heatOp);
    Vector<double> ones = heatSolver.solve(rho);
    Vector<double> scalarExtension = heatSolver.solve(rho_magnitudes);
    for (size_t i = 0; i < V; i++) {
        scalarExtension[i] /= ones[i];
    }

    // Vector extension
    SparseMatrix<std::complex<double>> CL = geometry->polygonVertexConnectionLaplacian;
    SparseMatrix<std::complex<double>> M_CL = M.cast<std::complex<double>>();
    SparseMatrix<std::complex<double>> vectorHeatOp = M_CL + shortTime * CL;
    Vector<std::complex<double>> vectorSoln = solvePositiveDefinite(vectorHeatOp, rho_vectors);
    VertexData<Vector2> vectorExtension(*mesh);
    for (Vertex v : mesh->vertices()) {
        size_t vIdx = geometry->vertexIndices[v];
        Vector2 vec = Vector2::fromComplex(vectorSoln[vIdx]);
        vectorExtension[v] = vec / vec.norm();
    }

    VertexData<Vector2> polySoln(*mesh);
    for (size_t i = 0; i < V; i++) {
        polySoln[i] = scalarExtension[i] * vectorExtension[i];
    }

    // Visualize solution.
    VertexData<Vector3> polySolnR3(*mesh);
    for (size_t i = 0; i < V; i++) {
        polySolnR3[i] = polySoln[i][0] * VBASISX[i] + polySoln[i][1] * VBASISY[i];
    }
    psMesh->addVertexVectorQuantity("vector sources", vectorSources);
    psMesh->addVertexVectorQuantity("VHM", polySolnR3);

    if (mesh->isTriangular()) {
        // Solve using standard operators.
        std::unique_ptr<VertexPositionGeometry> manifoldGeom;
        std::unique_ptr<ManifoldSurfaceMesh> manifoldMesh;
        manifoldMesh = mesh->toManifoldMesh();
        manifoldGeom = geometry->reinterpretTo(*manifoldMesh);

        manifoldGeom->requireVertexConnectionLaplacian();
        // std::cerr << geometry->polygonVertexConnectionLaplacian << std::endl;
        // std::cerr << manifoldGeom->vertexConnectionLaplacian << std::endl;
        // std::cerr << geometry->polygonVertexConnectionLaplacian - manifoldGeom->vertexConnectionLaplacian <<
        // std::endl;
        std::cerr << (geometry->polygonVertexConnectionLaplacian - manifoldGeom->vertexConnectionLaplacian).norm()
                  << "\t"
                  << (geometry->polygonVertexConnectionLaplacian - manifoldGeom->vertexConnectionLaplacian).norm() /
                         (manifoldGeom->vertexConnectionLaplacian).norm()
                  << std::endl;
        manifoldGeom->unrequireVertexConnectionLaplacian();


        double tCoef = (MAX_DIAGONAL_LENGTH * MAX_DIAGONAL_LENGTH) / (MEAN_EDGE_LENGTH * MEAN_EDGE_LENGTH);
        VectorHeatMethodSolver triSolver(*manifoldGeom, tCoef);
        // re-interpret sources on manifold mesh
        std::vector<std::tuple<Vertex, double>> sourceMagnitudesManifold(nSources);
        std::vector<std::tuple<Vertex, Vector2>> sourceVectorsManifold(nSources);
        for (int i = 0; i < nSources; i++) {
            size_t vIdx = geometry->vertexIndices[std::get<0>(sourceMagnitudes[i])];
            Vertex manifoldVertex = manifoldMesh->vertex(vIdx);
            sourceMagnitudesManifold[i] = std::make_tuple(manifoldVertex, std::get<1>(sourceMagnitudes[i]));
            sourceVectorsManifold[i] = std::make_tuple(manifoldVertex, std::get<1>(sourceVectors[i]));
        }
        VertexData<double> triScalarExtension = triSolver.extendScalar(sourceMagnitudesManifold);
        VertexData<Vector2> triVectorExtension = triSolver.transportTangentVectors(sourceVectorsManifold);
        VertexData<Vector2> triSoln(*manifoldMesh);
        for (Vertex v : manifoldMesh->vertices()) triSoln[v] = triVectorExtension[v] * triScalarExtension[v];

        psMesh->addVertexTangentVectorQuantity("VHM (tri)", triSoln, VBASISX, VBASISY);
        psMesh->addVertexTangentVectorQuantity("VHM (poly)", polySoln, VBASISX, VBASISY);
        psMesh->addVertexScalarQuantity("scalar ext. (tri)", triScalarExtension);
        psMesh->addVertexScalarQuantity("scalar ext. (poly)", scalarExtension);
        std::cerr << "|scalars - truth|: " << (triScalarExtension.toVector() - scalarExtension).norm() << "\t"
                  << (triScalarExtension.toVector() - scalarExtension).norm() / triScalarExtension.toVector().norm()
                  << std::endl;
        double vectorResidual = 0.;
        double solnResidual = 0.;
        for (Vertex v : mesh->vertices()) {
            size_t vIdx = geometry->vertexIndices[v];
            vectorResidual += (triVectorExtension[vIdx] - Vector2::fromComplex(vectorExtension[vIdx])).norm();
            solnResidual += (polySoln[v] - triSoln[vIdx]).norm();
        }
        std::cerr << "|vectors - truth|: " << vectorResidual << "\t" << vectorResidual / V << std::endl;
        std::cerr << "|soln - truth|: " << solnResidual << "\t" << solnResidual / V << std::endl;
    }

    geometry->unrequireVertexIndices();
    geometry->unrequirePolygonLaplacian();
    geometry->unrequirePolygonVertexConnectionLaplacian();
    geometry->unrequirePolygonVertexLumpedMassMatrix();
}

void myCallback() {

    if (ImGui::Button("Timing")) {
        timing();
    }
    if (ImGui::Button("Test mass lumping")) {
        testMassLumping();
    }
    if (mesh->isTriangular()) {
        if (ImGui::Button("Test on triangle meshes")) {
            testOnTriangleMeshes();
        }
    }
    if (ImGui::Button("Test DEC operators")) {
        testDECOperators();
    }
    if (ImGui::Button("Solve Poisson problems")) {
        solvePoissonProblems();
    }
    if (ImGui::Button("Solve geodesic distance")) {
        solveGeodesicDistance();
    }
    if (ImGui::Button("Solve Vector Heat Method")) {
        solveVectorHeatMethod();
    }
    if (ImGui::Button("Solve Signed Heat Method")) {
        VertexData<double> phi = signedHeatSolver->computeDistance(CURVES);
        psMesh->addVertexSignedDistanceQuantity("GSD", phi)->setEnabled(true);
    }
}

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("A program demonstrating use of various DEC operators for polygon meshes.");
    args::Positional<std::string> meshFilename(parser, "mesh", "A mesh file.");
    args::ValueFlag<std::string> inputFilename(parser, "input", "Input curve filepath", {"i", "input"});

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help&) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    polyscope::init();

    polyscope::state::userCallback = myCallback;

    // Load mesh
    if (meshFilename) {
        std::string MESH_FILEPATH = args::get(meshFilename);
        std::string MESHNAME = polyscope::guessNiceNameFromPath(MESH_FILEPATH);
        std::tie(mesh, geometry) = readSurfaceMesh(MESH_FILEPATH);
        psMesh = polyscope::registerSurfaceMesh(MESHNAME, geometry->vertexPositions, mesh->getFaceVertexList());
        // psMesh->setAllPermutations(polyscopePermutations(*mesh)); // not valid for polygon meshes

        MAX_DIAGONAL_LENGTH = 0.; // maximum length over all polygon diagonals
        MEAN_EDGE_LENGTH = 0.;
        for (Face f : mesh->faces()) {
            std::vector<Vertex> vertices;
            for (Vertex v : f.adjacentVertices()) {
                vertices.push_back(v);
            }
            size_t n = vertices.size();
            for (size_t i = 0; i < n; i++) {
                Vector3 pi = geometry->vertexPositions[vertices[i]];
                for (size_t j = i + 1; j < n; j++) {
                    Vector3 pj = geometry->vertexPositions[vertices[j]];
                    double length = (pi - pj).norm();
                    MAX_DIAGONAL_LENGTH = std::max(MAX_DIAGONAL_LENGTH, length);
                }
            }
        }
        geometry->requireEdgeLengths();
        for (Edge e : mesh->edges()) {
            MEAN_EDGE_LENGTH += geometry->edgeLengths[e];
        }
        MEAN_EDGE_LENGTH /= mesh->nEdges();
        geometry->unrequireEdgeLengths();

        VertexData<Vector3> vertexNormals(*mesh);
        if (mesh->isTriangular()) {
            geometry->requireVertexNormals();
            vertexNormals = geometry->vertexNormals;
            geometry->unrequireVertexNormals();
        } else {
            for (Vertex v : mesh->vertices()) {
                double totalArea = 0.;
                Vector3 vN = {0., 0., 0.};
                for (Face f : v.adjacentFaces()) {
                    // polygon vector area
                    Vector3 N = {0, 0, 0};
                    for (Halfedge he : f.adjacentHalfedges()) {
                        Vertex vA = he.vertex();
                        Vertex vB = he.next().vertex();
                        Vector3 pA = geometry->vertexPositions[vA];
                        Vector3 pB = geometry->vertexPositions[vB];
                        N += cross(pA, pB);
                    }
                    N *= 0.5;
                    vN += N;
                    totalArea += N.norm();
                }
                vN /= totalArea;
                vertexNormals[v] = vN;
            }
        }

        VBASISX.resize(mesh->nVertices());
        VBASISY.resize(mesh->nVertices());
        geometry->requireVertexIndices();
        for (Vertex v : mesh->vertices()) {
            size_t vIdx = geometry->vertexIndices[v];
            Vector3 xVec = geometry->halfedgeVector(v.halfedge());
            xVec /= xVec.norm();
            VBASISX[vIdx] = xVec;
            VBASISY[vIdx] = cross(vertexNormals[v], xVec);
        }
        geometry->requireVertexIndices();

        signedHeatSolver = std::unique_ptr<SignedHeatPolygon>(new SignedHeatPolygon(*geometry));
    }

    // Load source geometry.
    if (inputFilename) {
        std::string filename = args::get(inputFilename);
        CURVES = readCurves(filename);
        displayInput(CURVES);
    }

    polyscope::show();

    return EXIT_SUCCESS;
}