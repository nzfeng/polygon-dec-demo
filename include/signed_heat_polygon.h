#pragma once

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/surface_mesh.h"

#include "polyscope/curve_network.h"
#include "polyscope/surface_mesh.h"

#include <chrono>
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

using namespace geometrycentral;
using namespace geometrycentral::surface;

class SignedHeatPolygon {

  public:
    SignedHeatPolygon(EmbeddedGeometryInterface& geom, double tCoef = 1.0);

    VertexData<double> computeDistance(const std::vector<std::vector<Vertex>>& curves); // TODO: Add options

  private:
    SurfaceMesh& mesh;
    EmbeddedGeometryInterface& geom;

    double shortTime;
    double maxDiagonalLength;
    bool timeUpdated = false;

    // Solvers
    std::unique_ptr<PositiveDefiniteSolver<std::complex<double>>> vectorHeatSolver;
    std::unique_ptr<PositiveDefiniteSolver<double>> poissonSolver;
    SparseMatrix<double> laplaceMat;

    void ensureHaveVectorHeatSolver();
    void ensureHavePoissonSolver();

    void buildSignedCurveSource(const std::vector<Vertex>& curve, Vector<std::complex<double>>& X0) const;
};