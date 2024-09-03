#include "signed_heat_polygon.h"

SignedHeatPolygon::SignedHeatPolygon(EmbeddedGeometryInterface& geom_, double tCoef) : mesh(geom_.mesh), geom(geom_) {

    // Compute the maximum length over all polygon diagonals, as suggested by de Goes et al.
    maxDiagonalLength = 0.;
    for (Face f : mesh.faces()) {
        std::vector<Vertex> vertices;
        for (Vertex v : f.adjacentVertices()) {
            vertices.push_back(v);
        }
        size_t n = vertices.size();
        for (size_t i = 0; i < n; i++) {
            Vector3 pi = geom.vertexPositions[vertices[i]];
            for (size_t j = i + 1; j < n; j++) {
                Vector3 pj = geom.vertexPositions[vertices[j]];
                double length = (pi - pj).norm();
                maxDiagonalLength = std::max(maxDiagonalLength, length);
            }
        }
    }
    shortTime = tCoef * maxDiagonalLength;

    // geom.requireEdgeLengths();
    // double meanEdgeLength = 0.;
    // for (Edge e : mesh.edges()) {
    //     meanEdgeLength += geom.edgeLengths[e];
    // }
    // meanEdgeLength /= mesh.nEdges();
    // geom.unrequireEdgeLengths();
    // shortTime = tCoef * meanEdgeLength * meanEdgeLength;

    geom.requirePolygonLaplacian();
    laplaceMat = geom.polygonLaplacian;
    geom.unrequirePolygonLaplacian();
}

VertexData<double> SignedHeatPolygon::computeDistance(const std::vector<std::vector<Vertex>>& curves) {

    size_t V = mesh.nVertices();
    Vector<std::complex<double>> X0 = Vector<std::complex<double>>::Zero(V);
    geom.requireVertexIndices();
    geom.requireVertexTangentBasis();
    geom.requireVertexNormals();
    for (const auto& curve : curves) buildSignedCurveSource(curve, X0);
    geom.unrequireVertexNormals();

    ensureHaveVectorHeatSolver();
    Vector<std::complex<double>> Xt = vectorHeatSolver->solve(X0);
    // TODO: normals preservation

    // Average onto faces, and normalize.
    size_t F = mesh.nFaces();
    Vector<double> Y(3 * F);
    Eigen::MatrixXd X(F, 3); // TODO: remove
    geom.requireFaceIndices();
    for (Face f : mesh.faces()) {
        Vector3 Yf = {0, 0, 0};
        for (Vertex v : f.adjacentVertices()) {
            size_t vIdx = geom.vertexIndices[v];
            Yf += std::real(Xt[vIdx]) * geom.vertexTangentBasis[v][0];
            Yf += std::imag(Xt[vIdx]) * geom.vertexTangentBasis[v][1];
        }
        Yf /= Yf.norm();
        size_t fIdx = geom.faceIndices[f];
        for (int j = 0; j < 3; j++) {
            Y[3 * fIdx + j] = Yf[j];
            X(fIdx, j) = Yf[j]; // TODO: remove
        }
    }
    geom.unrequireFaceIndices();
    geom.unrequireVertexTangentBasis();

    VertexData<Vector3> vBasisX(mesh);
    VertexData<Vector3> vBasisY(mesh);
    for (Vertex v : mesh.vertices()) {
        vBasisX[v] = geom.vertexTangentBasis[v][0];
        vBasisY[v] = geom.vertexTangentBasis[v][1];
    }
    polyscope::getSurfaceMesh("polygon-bear")->addVertexTangentVectorQuantity("Xt", Xt, vBasisX, vBasisY);
    polyscope::getSurfaceMesh("polygon-bear")->addVertexTangentVectorQuantity("X0", X0, vBasisX, vBasisY);
    polyscope::getSurfaceMesh("polygon-bear")->addFaceVectorQuantity("Y", X); // TODO: remove

    geom.requirePolygonDivergenceMatrix();
    Vector<double> divYt = geom.polygonDivergenceMatrix * Y;
    geom.unrequirePolygonDivergenceMatrix();

    Vector<double> phi;
    Vector<bool> setAMembership = Vector<bool>::Ones(V);
    for (const auto& curve : curves) {
        for (const Vertex& v : curve) {
            setAMembership[geom.vertexIndices[v]] = false;
        }
    }
    int nB = V - setAMembership.cast<int>().sum();
    Vector<double> bcVals = Vector<double>::Zero(nB);
    BlockDecompositionResult<double> decomp = blockDecomposeSquare(laplaceMat, setAMembership, true);
    Vector<double> rhsValsA, rhsValsB;
    decomposeVector(decomp, divYt, rhsValsA, rhsValsB);
    Vector<double> combinedRHS = rhsValsA;
    Vector<double> Aresult = solvePositiveDefinite(decomp.AA, combinedRHS);
    phi = reassembleVector(decomp, Aresult, bcVals);

    geom.unrequireVertexIndices();

    // TODO: level set constraints

    return VertexData<double>(mesh, phi);
}

void SignedHeatPolygon::buildSignedCurveSource(const std::vector<Vertex>& curve,
                                               Vector<std::complex<double>>& X0) const {

    // Encode curve input by expressing curve normals in vertex tangent spaces.
    size_t nNodes = curve.size();
    for (size_t i = 0; i < nNodes - 1; i++) {
        Vertex vA = curve[i];
        Vertex vB = curve[i + 1];
        size_t vIdxA = geom.vertexIndices[vA];
        size_t vIdxB = geom.vertexIndices[vB];
        Vector3 pA = geom.vertexPositions[vA];
        Vector3 pB = geom.vertexPositions[vB];
        Vector3 tangent = pB - pA;
        Vector3 vnA = geom.vertexNormals[vA];
        Vector3 vnB = geom.vertexNormals[vB];
        Vector3 xAxisA = geom.vertexTangentBasis[vA][0];
        Vector3 yAxisA = geom.vertexTangentBasis[vA][1];
        Vector3 xAxisB = geom.vertexTangentBasis[vB][0];
        Vector3 yAxisB = geom.vertexTangentBasis[vB][1];
        Vector3 tangentA = dot(xAxisA, tangent) * xAxisA + dot(yAxisA, tangent) * yAxisA;
        Vector3 tangentB = dot(xAxisB, tangent) * xAxisB + dot(yAxisB, tangent) * yAxisB;
        Vector3 normalA = cross(vnA, tangentA);
        Vector3 normalB = cross(vnB, tangentB);
        X0[vIdxA] += std::complex<double>(dot(geom.vertexTangentBasis[vA][0], normalA),
                                          dot(geom.vertexTangentBasis[vA][1], normalA));
        X0[vIdxB] += std::complex<double>(dot(geom.vertexTangentBasis[vB][0], normalB),
                                          dot(geom.vertexTangentBasis[vB][1], normalB));
    }
}

void SignedHeatPolygon::ensureHaveVectorHeatSolver() {

    if (vectorHeatSolver != nullptr && !timeUpdated) return;
    timeUpdated = false;

    geom.requirePolygonVertexConnectionLaplacian();
    geom.requirePolygonVertexLumpedMassMatrix();

    SparseMatrix<std::complex<double>>& Lconn = geom.polygonVertexConnectionLaplacian;
    SparseMatrix<double>& massMat = geom.polygonVertexLumpedMassMatrix;
    SparseMatrix<std::complex<double>> vectorOp = massMat.cast<std::complex<double>>() + shortTime * Lconn;
    vectorHeatSolver.reset(new PositiveDefiniteSolver<std::complex<double>>(vectorOp));

    geom.unrequirePolygonVertexConnectionLaplacian();
    geom.unrequirePolygonVertexLumpedMassMatrix();
}

void SignedHeatPolygon::ensureHavePoissonSolver() {

    if (poissonSolver != nullptr) return;

    poissonSolver.reset(new PositiveDefiniteSolver<double>(laplaceMat));
}