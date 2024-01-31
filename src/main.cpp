#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

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
    // TODO

    if (mesh->isTriangular()) {
        // TODO: Solve using standard cotan operator.
    }
}

/* Check that gradient, divergence, Laplacian can be assembled as expected from DEC operators. */
void testDECOperators() {
    // TODO
}

/* Implement the scalar heat method -- a great way to test the gradient, divergence, and Laplacian. */
void solveGeodesicDistance() {
    // TODO
}

/* Test that the lumped mass matrices correspond with their unlumped versions. */
void testMassLumping() {}

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

    polyscope::show();

    return EXIT_SUCCESS;
}