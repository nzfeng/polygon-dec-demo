#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/utilities.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

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

/* Solve Poisson problems using both sets of polygon operators. */
void solvePoissonProblems() {
    std::cerr << "solvePoissonProblems()..." << std::endl;
    // Randomly choose some vertex sources.
    size_t V = mesh->nVertices();
    int nSources = randomInt(1, 5);
    Vector<double> rho = Vector<double>::Zero(V);
    for (int i = 0; i < nSources; i++) {
        size_t vIdx = randomIndex(V);
        double weight = randomInt(1, 5);
        weight *= (unitRand() < 0.5) ? 1. : -1.;
        rho[vIdx] += weight;
    }

    geometry->requireVirtualRefinementLaplacian();
    geometry->requireVirtualRefinementVertexLumpedMassMatrix();
    geometry->requireVirtualElementLaplacian();
    geometry->requireVirtualElementVertexLumpedMassMatrix();

    SparseMatrix<double> L_VR = geometry->virtualRefinementLaplacian;
    SparseMatrix<double> L_VEM = geometry->virtualElementLaplacian;
    SparseMatrix<double> M_VR = geometry->virtualRefinementVertexLumpedMassMatrix;
    SparseMatrix<doulbe> M_VEM = geometry->virtualElementVertexLumpedMassMatrix;
    double totalRho_VR = (M_VR * rho).sum();
    double totalRho_VEM = (M_VEM * rho).sum();
    double totalArea_VR = 0.;
    double totalArea_VEM = 0.;
    Vector<double> rhoBar_VR = Vector<double>::Ones(V) * (totalRho_VR / totalArea_VR);
    Vector<double> rhs_VR = M_VR * (rhoBar_VR - rho_VR);
    Vector<double> rhoBar_VEM = Vector<double>::Ones(V) * (totalRho_VEM / totalArea_VEM);
    Vector<double> rhs_VEM = M_VEM * (rhoBar_VEM - rho_VEM);
    Vector<double> VRSoln = solvePositiveDefinite(L_VR, rhs_VR);
    Vector<double> VEMSoln = solvePositiveDefinite(L_VEM, rhs_VEM);
    std::cerr << "|VRSoln - VEMSoln| / |VRSoln|: " << (VRSoln - VEMSoln).norm() / VRSoln.norm() << std::endl;
    std::cerr << "|VRSoln - VEMSoln| / |VEMSoln|: " << (VRSoln - VEMSoln).norm() / VEMSoln.norm() << std::endl;

    if (mesh->isTriangular()) {
        // Solve using standard cotan operator.
        geometry->requireVertexLumpedMassMatrix();
        geometry->requireCotanLaplacian();
        geometry->requireFaceAreas();

        SparseMatrix<double> C = geometry->cotanLaplacian;
        SparseMatrix<double> M = geometry->vertexLumpedMassMatrix;
        double totalArea = 0.;
        for (Face f : mesh->faces()) totalArea += geometry->faceAreas[f];
        double totalRho = (M * rho).sum();
        Vector<double> rhoBar = Vector<double>::Ones(V) * (totalRho / totalArea);
        Vector<double> rhs = M * (rhoBar - rho);
        Vector<double> triSoln = solvePositiveDefinite(C, rhs);
        psMesh->addVertexSignedDistanceQuantity("tri", triSoln);
        std::cerr << "[poisson] |VR - truth|: " << (triSoln - VRSoln).norm() << std::endl;
        std::cerr << "[poisson] |VEM - truth|: " << (triSoln - VEMSoln).norm() << std::endl;

        geometry->unrequireCotanLaplacian();
        geometry->unrequireVertexLumpedMassMatrix();
        geometry->unrequireFaceAreas();
    }

    geometry->unrequireVirtualRefinementLaplacian();
    geometry->unrequireVirtualRefinementVertexLumpedMassMatrix();
    geometry->unrequireVirtualElementLaplacian();
    geometry->unrequireVirtualElementVertexLumpedMassMatrix();
    std::cerr << "\tDone testing." << std::endl;
}

/* Check that gradient, divergence, Laplacian can be assembled as expected from DEC operators. */
void testDECOperators() {
    std::cerr << "testDECOperators()..." << std::endl;
    geometry->requireVirtualRefinementDECOperators();
    geometry->requireVirtualRefinementLaplacian();
    geometry->requireVirtualElementDECOperators();
    geometry->requireVirtualElementLaplacian();

    // TODO
    double epsilon = 1e-8;
    SparseMatrix<double>& L_VR = geometry->virtualRefinementLaplacian;
    SparseMatrix<double>& L_VEM = geometry->virtualElementLaplacian;
    SparseMatrix<double>& h0_VR = geometry->virtualRefinementHodge0;
    SparseMatrix<double>& h0_VEM = geometry->virtualElementHodge0;
    SparseMatrix<double>& h0Inv_VR = geometry->virtualRefinementHodge0Inverse;
    SparseMatrix<double>& h0Inv_VEM = geometry->virtualElementHodge0Inverse;
    SparseMatrix<double>& h1_VR = geometry->virtualRefinementHodge1;
    SparseMatrix<double>& h1_VEM = geometry->virtualElementHodge1;
    SparseMatrix<double>& h1Inv_VR = geometry->virtualRefinementHodge1Inverse;
    SparseMatrix<double>& h1Inv_VEM = geometry->virtualElementHodge1Inverse;
    SparseMatrix<double>& h2_VR = geometry->virtualRefinementHodge2;
    SparseMatrix<double>& h2_VEM = geometry->virtualElementHodge2;
    SparseMatrix<double>& h2Inv_VR = geometry->virtualRefinementHodge2Inverse;
    SparseMatrix<double>& h2Inv_VEM = geometry->virtualElementHodge2Inverse;
    SparseMatrix<double>& d0_VR = geometry->virtualRefinementD0;
    SparseMatrix<double>& d0_VEM = geometry->virtualElementD0;
    SparseMatrix<double>& d1_VR = geometry->virtualRefinementD1;
    SparseMatrix<double>& d1_VEM = geometry->virtualElementD1;
    assert((L_VR - d0_VR.transpose() * h1_VR * d0_VR).norm() < epsilon);
    assert((L_VEM - d0_VEM.transpose() * h1_VEM * d0_VEM).norm() < epsilon);

    geometry->unrequireVirtualRefinementDECOperators();
    geometry->unrequireVirtualRefinementLaplacian();
    geometry->unrequireVirtualElementDECOperators();
    geometry->unrequireVirtualElementLaplacian();
    std::cerr << "\tDone testing." << std::endl;
}

/* Implement the scalar heat method -- good way to test the gradient, divergence, and Laplacian. */
void solveGeodesicDistance() {
    std::cerr << "solveGeodesicDistance()..." << std::endl;
    // Randomly choose some vertex sources.
    size_t V = mesh->nVertices();
    int nSources = randomInt(1, 5);
    std::vector<Vertex> sources(nSources);
    for (int i = 0; i < nSources; i++) {
        sources[i] = mesh->vertex(randomIndex(V));
    }

    // TODO
    Vector<double> VRDistances, VEMDistances;

    if (mesh->isTriangular()) {
        // Solve using standard operators.
        HeatMethodDistanceSolver triSolver(*geometry);
        VertexData<double> triDistances = triSolver.computeDistance(sources);
        psMesh->addVertexSignedDistanceQuantity("tri", triDistances);
        std::cerr << "[distance] |VR - truth|: " << (triDistances - VRDistances).norm() << std::endl;
        std::cerr << "[distance] |VEM - truth|: " << (triDistances - VEMDistances).norm() << std::endl;
    }
    std::cerr << "\tDone testing." << std::endl;
}

/* Test that the lumped mass matrices correspond with their unlumped versions. */
void testMassLumping() {
    std::cerr << "testMassLumping()..." << std::endl;
    geometry->requireVirtualRefinementVertexLumpedMassMatrix();
    geometry->requireVirtualRefinementVertexGalerkinMassMatrix();
    geometry->requireVirtualElementVertexLumpedMassMatrix();
    geometry->requireVirtualElementVertexGalerkinMassMatrix();

    double epsilon = 1e-8;
    SparseMatrix<double> M_VR_T = geometry->virtualRefinementVertexGalerkinMassMatrix.transpose();
    SparseMatrix<double> M_VEM_T = geometry->virtualElementVertexGalerkinMassMatrix.transpose();
    SparseMatrix<double> M_VR = geometry->virtualRefinementVertexLumpedMassMatrix;
    SparseMatrix<double> M_VEM = geometry->virtualElementVertexLumpedMassMatrix;
    for (size_t i = 0; i < mesh->nVertices(); i++) {
        double rowSumVR = 0.;
        double rowSumVEM = 0.;
        for (SparseMatrix<double>::InnerIterator it(M_VR_T, i); it; ++it) {
            rowSumVR += it.value();
        }
        for (SparseMatrix<double>::InnerIterator it(M_VEM_T, i); it; ++it) {
            rowSumVEM += it.value();
        }
        assert(std::abs(rowSumVR - M_VR.coeffRef(i, i)) < epsilon);
        assert(std::abs(rowSumVEM - M_VEM.coeffRef(i, i)) < epsilon);
    }

    geometry->unrequireVirtualRefinementVertexLumpedMassMatrix();
    geometry->unrequireVirtualRefinementVertexGalerkinMassMatrix();
    geometry->unrequireVirtualElementVertexLumpedMassMatrix();
    geometry->unrequireVirtualElementVertexGalerkinMassMatrix();
    std::cerr << "\tDone testing." << std::endl;
}

void solveVectorHeatMethod() {

    // Randomly choose some vertex sources with random magnitudes.
    size_t V = mesh->nVertices();
    int nSources = randomInt(1, 5);
    std::vector<std::tuple<Vertex, double>> sourceMagnitudes(nSources);
    std::vector<std::tuple<Vertex, Vector2>> sourceVectors(nSources);
    for (int i = 0; i < nSources; i++) {
        points[i] = randomReal(1., 5.);
        Vector2 vec = {randomReal(-1., 1.), randomReal(-1., 1.)};
        sourceVectors[i] = vec / vec.norm();
    }

    // TODO

    if (mesh->isTriangular()) {
        // Solve using standard operators.
        VectorHeatMethodSolver triSolver(*geometry);
        VertexData<double> triScalarExtension = triSolver.extendScalar(sourceMagnitudes);
        VertexData<Vector2> = triSolver.transportTangentVectors(sourceVectors);
        // TODO
    }
}

void myCallback() {}

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("A program demonstrating use of various DEC operators for polygon meshes.");
    args::Positional<std::string> meshFilename(parser, "mesh", "A mesh file.");

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
        psMesh->setAllPermutations(polyscopePermutations(*mesh));
    }

    solvePoissonProblems();
    testDECOperators();
    solveGeodesicDistance();
    testMassLumping();

    polyscope::show();

    return EXIT_SUCCESS;
}