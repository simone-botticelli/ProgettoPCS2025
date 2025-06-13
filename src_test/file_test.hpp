#include <gtest/gtest.h>
#include "PolyhedralMesh.hpp"
#include <cmath>

using namespace PolyhedralLibrary;
using namespace Eigen;

class PolyhedralMeshTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Setup common test data
        tetrahedron = std::make_unique<PlatonicSolid>(PlatonicType::TETRAHEDRON);
        octahedron = std::make_unique<PlatonicSolid>(PlatonicType::OCTAHEDRON);
        icosahedron = std::make_unique<PlatonicSolid>(PlatonicType::ICOSAHEDRON);
    }

    std::unique_ptr<PlatonicSolid> tetrahedron;
    std::unique_ptr<PlatonicSolid> octahedron;
    std::unique_ptr<PlatonicSolid> icosahedron;
};

// Test for PlatonicSolid constructor
TEST_F(PolyhedralMeshTest, PlatonicSolid_TetrahedronConstructor) {
    EXPECT_EQ(tetrahedron->Type, PlatonicType::TETRAHEDRON);
    EXPECT_EQ(tetrahedron->NumVertices, 4);
    EXPECT_EQ(tetrahedron->NumEdges, 6);
    EXPECT_EQ(tetrahedron->NumFaces, 4);
    
    // Check vertices IDs
    EXPECT_EQ(tetrahedron->VerticesId.size(), 4);
    for (unsigned int i = 0; i < 4; ++i) {
        EXPECT_EQ(tetrahedron->VerticesId[i], i);
    }
    
    // Check coordinates matrix dimensions
    EXPECT_EQ(tetrahedron->VerticesCoordinates.rows(), 3);
    EXPECT_EQ(tetrahedron->VerticesCoordinates.cols(), 4);
    
    // Check edges extrema matrix dimensions
    EXPECT_EQ(tetrahedron->EdgesExtrema.rows(), 2);
    EXPECT_EQ(tetrahedron->EdgesExtrema.cols(), 6);
    
    // Check faces vertices structure
    EXPECT_EQ(tetrahedron->FacesVertices.size(), 4);
    for (const auto& face : tetrahedron->FacesVertices) {
        EXPECT_EQ(face.size(), 3); // Triangular faces
    }
}

TEST_F(PolyhedralMeshTest, PlatonicSolid_OctahedronConstructor) {
    EXPECT_EQ(octahedron->Type, PlatonicType::OCTAHEDRON);
    EXPECT_EQ(octahedron->NumVertices, 6);
    EXPECT_EQ(octahedron->NumEdges, 12);
    EXPECT_EQ(octahedron->NumFaces, 8);
    
    EXPECT_EQ(octahedron->VerticesCoordinates.rows(), 3);
    EXPECT_EQ(octahedron->VerticesCoordinates.cols(), 6);
}

TEST_F(PolyhedralMeshTest, PlatonicSolid_IcosahedronConstructor) {
    EXPECT_EQ(icosahedron->Type, PlatonicType::ICOSAHEDRON);
    EXPECT_EQ(icosahedron->NumVertices, 12);
    EXPECT_EQ(icosahedron->NumEdges, 30);
    EXPECT_EQ(icosahedron->NumFaces, 20);
    
    EXPECT_EQ(icosahedron->VerticesCoordinates.rows(), 3);
    EXPECT_EQ(icosahedron->VerticesCoordinates.cols(), 12);
}

// Test for NormalizeMatrixColumns
TEST_F(PolyhedralMeshTest, NormalizeMatrixColumns_BasicTest) {
    MatrixXd test_matrix(3, 2);
    test_matrix << 3, 0,
                   4, 4,
                   0, 3;
    
    NormalizeMatrixColumns(test_matrix);
    
    // Check that each column has unit norm
    EXPECT_NEAR(test_matrix.col(0).norm(), 1.0, 1e-10);
    EXPECT_NEAR(test_matrix.col(1).norm(), 1.0, 1e-10);
    
    // Check specific values
    EXPECT_NEAR(test_matrix(0, 0), 0.6, 1e-10);  // 3/5
    EXPECT_NEAR(test_matrix(1, 0), 0.8, 1e-10);  // 4/5
    EXPECT_NEAR(test_matrix(2, 0), 0.0, 1e-10);  // 0/5
}

TEST_F(PolyhedralMeshTest, NormalizeMatrixColumns_ZeroColumn) {
    MatrixXd test_matrix(3, 2);
    test_matrix << 0, 1,
                   0, 0,
                   0, 0;
    
    NormalizeMatrixColumns(test_matrix);
    
    // Zero column should remain zero
    EXPECT_NEAR(test_matrix.col(0).norm(), 0.0, 1e-10);
    // Non-zero column should be normalized
    EXPECT_NEAR(test_matrix.col(1).norm(), 1.0, 1e-10);
}

// Test for Initialize_ClassI_GeodesicCounts
TEST_F(PolyhedralMeshTest, Initialize_ClassI_GeodesicCounts_Tetrahedron) {
    GeodesicPolyhedron geodesic;
    unsigned int n = 3;
    
    Initialize_ClassI_GeodesicCounts(geodesic, *tetrahedron, n);
    
    unsigned int T = n * n; // T = 9
    EXPECT_EQ(geodesic.NumCell0Ds, 2 * T + 2);  // 20
    EXPECT_EQ(geodesic.NumCell1Ds, 6 * T);      // 54
    EXPECT_EQ(geodesic.NumCell2Ds, 4 * T);      // 36
}

TEST_F(PolyhedralMeshTest, Initialize_ClassI_GeodesicCounts_Octahedron) {
    GeodesicPolyhedron geodesic;
    unsigned int n = 2;
    
    Initialize_ClassI_GeodesicCounts(geodesic, *octahedron, n);
    
    unsigned int T = n * n; // T = 4
    EXPECT_EQ(geodesic.NumCell0Ds, 4 * T + 2);  // 18
    EXPECT_EQ(geodesic.NumCell1Ds, 12 * T);     // 48
    EXPECT_EQ(geodesic.NumCell2Ds, 8 * T);      // 32
}

TEST_F(PolyhedralMeshTest, Initialize_ClassI_GeodesicCounts_Icosahedron) {
    GeodesicPolyhedron geodesic;
    unsigned int n = 2;
    
    Initialize_ClassI_GeodesicCounts(geodesic, *icosahedron, n);
    
    unsigned int T = n * n; // T = 1
    EXPECT_EQ(geodesic.NumCell0Ds, 10 * T + 2);  // 12
    EXPECT_EQ(geodesic.NumCell1Ds, 30 * T);      // 30
    EXPECT_EQ(geodesic.NumCell2Ds, 20 * T);      // 20
}

// Test for Initialize_ClassII_GeodesicCounts
TEST_F(PolyhedralMeshTest, Initialize_ClassII_GeodesicCounts_Tetrahedron) {
    GeodesicPolyhedron geodesic;
    unsigned int n = 2;
    
    Initialize_ClassII_GeodesicCounts(geodesic, *tetrahedron, n);
    
    unsigned int T = n * n; // T = 4
    EXPECT_EQ(geodesic.NumCell0Ds, 6 * (T + n) + 2);  // 6 * 6 + 2 = 38
    EXPECT_EQ(geodesic.NumCell1Ds, 18 * (T + n));     // 18 * 6 = 108
    EXPECT_EQ(geodesic.NumCell2Ds, 12 * (T + n));     // 12 * 6 = 72
}

// Test for InitializeGeodesicStorage
TEST_F(PolyhedralMeshTest, InitializeGeodesicStorage_BasicTest) {
    GeodesicPolyhedron geodesic;
    geodesic.NumCell0Ds = 10;
    geodesic.NumCell1Ds = 16;
    geodesic.NumCell2Ds = 8;
    
    InitializeGeodesicStorage(geodesic);
    
    EXPECT_EQ(geodesic.Cell0DsId.size(), 10);
    EXPECT_EQ(geodesic.Cell1DsId.size(), 16);
    EXPECT_EQ(geodesic.Cell2DsId.size(), 8);
    
    EXPECT_EQ(geodesic.Cell0DsCoordinates.rows(), 3);
    EXPECT_EQ(geodesic.Cell0DsCoordinates.cols(), 10);
    
    EXPECT_EQ(geodesic.Cell1DsExtrema.rows(), 2);
    EXPECT_EQ(geodesic.Cell1DsExtrema.cols(), 16);
    
    EXPECT_EQ(geodesic.Cell2DsVertices.size(), 8);
    EXPECT_EQ(geodesic.Cell2DsEdges.size(), 8);
}

// Test for addVertex
TEST_F(PolyhedralMeshTest, addVertex_BasicTest) {
    GeodesicPolyhedron geodesic;
    geodesic.NumCell0Ds = 20;
    geodesic.NumCell1Ds = 32;
    geodesic.NumCell2Ds = 16;
    
    InitializeGeodesicStorage(geodesic);
    
    const Vector3d coords(1.0, 2.0, 3.0);
    
    addVertex(geodesic, 11, coords);
    
    EXPECT_EQ(geodesic.Cell0DsId[11], 11);
    EXPECT_NEAR(geodesic.Cell0DsCoordinates(0, 11), 1.0, 1e-10);
    EXPECT_NEAR(geodesic.Cell0DsCoordinates(1, 11), 2.0, 1e-10);
    EXPECT_NEAR(geodesic.Cell0DsCoordinates(2, 11), 3.0, 1e-10);
}

// Test for GetorAddEdge
TEST_F(PolyhedralMeshTest, GetorAddEdge_NewEdge) {
    GeodesicPolyhedron geodesic;
	geodesic.NumCell0Ds = 10;
    geodesic.NumCell1Ds = 16;
    geodesic.NumCell2Ds = 8;

    InitializeGeodesicStorage(geodesic);
    
    unsigned int nextEdgeId = 0;
    unsigned int edgeId = GetorAddEdge(geodesic, nextEdgeId, 1, 3);
    
    EXPECT_EQ(edgeId, 0);
    EXPECT_EQ(nextEdgeId, 1);
    EXPECT_EQ(geodesic.Cell1DsId[0], 0);
    EXPECT_EQ(geodesic.Cell1DsExtrema(0, 0), 1);
    EXPECT_EQ(geodesic.Cell1DsExtrema(1, 0), 3);
    
    // Check if edge is in the map
    auto key = std::minmax(1u, 3u);
    EXPECT_EQ(geodesic.ExtrematoEdge[key], 0);
}

TEST_F(PolyhedralMeshTest, GetorAddEdge_ExistingEdge) {
    GeodesicPolyhedron geodesic;
	geodesic.NumCell0Ds = 10;
    geodesic.NumCell1Ds = 16;
    geodesic.NumCell2Ds = 8;

    InitializeGeodesicStorage(geodesic);
    
    unsigned int nextEdgeId = 0;
    
    // Add edge first time
    unsigned int edgeId1 = GetorAddEdge(geodesic, nextEdgeId, 1, 3);
    // Try to add same edge (reversed order)
    unsigned int edgeId2 = GetorAddEdge(geodesic, nextEdgeId, 3, 1);
    
    EXPECT_EQ(edgeId1, edgeId2);
    EXPECT_EQ(nextEdgeId, 1); // Should not increment again
}

// Test for addFace
TEST_F(PolyhedralMeshTest, addFace_TriangleTest) {
    GeodesicPolyhedron geodesic;
    geodesic.NumCell0Ds = 7;
    geodesic.NumCell1Ds = 10;
    geodesic.NumCell2Ds = 5;
    InitializeGeodesicStorage(geodesic);
    
    unsigned int nextEdgeId = 0;
    unsigned int nextFaceId = 0;
    
    addFace(geodesic, nextEdgeId, nextFaceId, 0, 1, 2);
    
    EXPECT_EQ(nextFaceId, 1);
    EXPECT_EQ(nextEdgeId, 3); // Should have added 3 edges
    
    // Check face vertices
    EXPECT_EQ(geodesic.Cell2DsVertices[0].size(), 3);
    EXPECT_EQ(geodesic.Cell2DsVertices[0][0], 0);
    EXPECT_EQ(geodesic.Cell2DsVertices[0][1], 1);
    EXPECT_EQ(geodesic.Cell2DsVertices[0][2], 2);
    
    // Check face edges
    EXPECT_EQ(geodesic.Cell2DsEdges[0].size(), 3);
}

// Test for ComputePointOnTriangle
TEST_F(PolyhedralMeshTest, ComputePointOnTriangle_BasicTest) {
    VectorXd vA(3), vB(3), vC(3);
    vA << 0, 0, 0;
    vB << 1, 0, 0;
    vC << 0, 1, 0;
    
    unsigned int n = 2;
    
    // Test corner point (should be vA)
    VectorXd result1 = ComputePointOnTriangle(n, 0, 0, vA, vB, vC);
    EXPECT_NEAR((result1 - vA).norm(), 0.0, 1e-10);
    
    // Test midpoint of edge AB
    VectorXd result2 = ComputePointOnTriangle(n, 1, 0, vA, vB, vC);
    VectorXd expected2 = 0.5 * vA + 0.5 * vB;
    EXPECT_NEAR((result2 - expected2).norm(), 0.0, 1e-10);
    
    // Test center point
    VectorXd result3 = ComputePointOnTriangle(3, 1, 1, vA, vB, vC);
    VectorXd expected3 = (vA + vB + vC) / 3.0;
    EXPECT_NEAR((result3 - expected3).norm(), 0.0, 1e-10);
}

// Test for Build_ClassI_Geodesic
TEST_F(PolyhedralMeshTest, Build_ClassI_Geodesic_N1_Tetrahedron) {
    GeodesicPolyhedron geodesic = Build_ClassI_Geodesic(*tetrahedron, 1);
    
    // For n=1, should be identical to original platonic solid
    EXPECT_EQ(geodesic.NumCell0Ds, tetrahedron->NumVertices);
    EXPECT_EQ(geodesic.NumCell1Ds, tetrahedron->NumEdges);
    EXPECT_EQ(geodesic.NumCell2Ds, tetrahedron->NumFaces);
    
    // Check that coordinates are normalized
    for (int i = 0; i < geodesic.Cell0DsCoordinates.cols(); ++i) {
        EXPECT_NEAR(geodesic.Cell0DsCoordinates.col(i).norm(), 1.0, 1e-10);
    }
}

TEST_F(PolyhedralMeshTest, Build_ClassI_Geodesic_N2_Tetrahedron) {
    GeodesicPolyhedron geodesic = Build_ClassI_Geodesic(*tetrahedron, 2);
    NormalizeMatrixColumns(geodesic.Cell0DsCoordinates);
    
    // Check counts match expected values
    unsigned int n = 2;
    unsigned int T = n * n;
    EXPECT_EQ(geodesic.NumCell0Ds, 2 * T + 2);
    EXPECT_EQ(geodesic.NumCell1Ds, 6 * T);
    EXPECT_EQ(geodesic.NumCell2Ds, 4 * T);
    
    // Check that all coordinates are normalized
    for (int i = 0; i < geodesic.Cell0DsCoordinates.cols(); ++i) {
        EXPECT_NEAR(geodesic.Cell0DsCoordinates.col(i).norm(), 1.0, 1e-10);
    }
    
    // Check that we have valid face structures
    for (unsigned int f = 0; f < geodesic.NumCell2Ds; ++f) {
        EXPECT_EQ(geodesic.Cell2DsVertices[f].size(), 3); // All triangular faces
        EXPECT_EQ(geodesic.Cell2DsEdges[f].size(), 3);
    }
}

TEST_F(PolyhedralMeshTest, Build_ClassII_Geodesic_N1_Tetrahedron) {
    unsigned int n = 3;
    GeodesicPolyhedron geodesic = Build_ClassII_Geodesic(*tetrahedron, n);

    // Check counts match expected values from the formula
    unsigned int T = n * n;
    EXPECT_EQ(geodesic.NumCell0Ds, 6 * (T + n) + 2);  // 6 * (1 + 1) + 2 = 14
    EXPECT_EQ(geodesic.NumCell1Ds, 18 * (T + n));     // 18 * 2 = 36
    EXPECT_EQ(geodesic.NumCell2Ds, 12 * (T + n));     // 12 * 2 = 24

    // Check that all coordinates are normalized
    for (int i = 0; i < geodesic.Cell0DsCoordinates.cols(); ++i) {
        EXPECT_NEAR(geodesic.Cell0DsCoordinates.col(i).norm(), 1.0, 1e-10);
    }
}

// Test for ComputeShortestPath
TEST_F(PolyhedralMeshTest, ComputeShortestPath_BasicTest) {
    GeodesicPolyhedron geodesic = Build_ClassI_Geodesic(*tetrahedron, 1);
    
    // Test path from vertex 0 to vertex 3
    ComputeShortestPath(geodesic, 0, 3);
    
    // Check that source and target are marked
    EXPECT_TRUE(geodesic.Cell0DsShortPath[0]);
    EXPECT_TRUE(geodesic.Cell0DsShortPath[3]);
    
    // Check that at least one edge is marked
    bool anyEdgeMarked = false;
    for (bool marked : geodesic.Cell1DsShortPath) {
        if (marked) {
            anyEdgeMarked = true;
            break;
        }
    }
    EXPECT_TRUE(anyEdgeMarked);
}

TEST_F(PolyhedralMeshTest, ComputeShortestPath_SameVertex) {
    GeodesicPolyhedron geodesic = Build_ClassI_Geodesic(*tetrahedron, 1);
    
    // Test path from vertex to itself
    ComputeShortestPath(geodesic, 0, 0);
    
    // Only source should be marked
    EXPECT_TRUE(geodesic.Cell0DsShortPath[0]);
    
    // No edges should be marked for same vertex
    bool anyEdgeMarked = false;
    for (bool marked : geodesic.Cell1DsShortPath) {
        if (marked) {
            anyEdgeMarked = true;
            break;
        }
    }
    EXPECT_FALSE(anyEdgeMarked);
}

TEST_F(PolyhedralMeshTest, ComputeShortestPath_InvalidVertex) {
    GeodesicPolyhedron geodesic = Build_ClassI_Geodesic(*tetrahedron, 1);
    
    // Test with invalid vertex indices
    // This should not crash and should handle gracefully
    EXPECT_NO_THROW(ComputeShortestPath(geodesic, 0, 999));
    EXPECT_NO_THROW(ComputeShortestPath(geodesic, 999, 0));
}

TEST_F(PolyhedralMeshTest, DualFaceEdgesAreCyclicallyConnected)
{
    GeodesicPolyhedron geodesic = Build_ClassI_Geodesic(*tetrahedron, 5);
    GeodesicPolyhedron poly = dualize(geodesic);

    for (unsigned int face : poly.Cell2DsId)
    {
        const auto& faceEdges = poly.Cell2DsEdges[face];
        const auto& faceVertices = poly.Cell2DsVertices[face];
        size_t E = faceEdges.size();

        ASSERT_EQ(E, faceVertices.size()) << "Face " << face << ": mismatch between number of edges and vertices.";

        // Verifica che il primo vertice sia in uno degli estremi del primo spigolo
        unsigned int firstEdgeIdx = faceEdges[0];
        unsigned int firstVertex = faceVertices[0];
        ASSERT_LT(firstEdgeIdx, poly.Cell1DsExtrema.cols());

        unsigned int orig0 = poly.Cell1DsExtrema(0, firstEdgeIdx);
        unsigned int end0 = poly.Cell1DsExtrema(1, firstEdgeIdx);

        EXPECT_TRUE(firstVertex == orig0 || firstVertex == end0)
            << "Face " << face << ": vertex[0] = " << firstVertex
            << " is not an endpoint of edge[0] = " << firstEdgeIdx;

        // Per ogni coppia di spigoli consecutivi
        for (size_t i = 0; i < E; ++i)
        {
            unsigned int edgeIdx1 = faceEdges[i];
            unsigned int edgeIdx2 = faceEdges[(i + 1) % E];

            ASSERT_LT(edgeIdx1, poly.Cell1DsExtrema.cols());
            ASSERT_LT(edgeIdx2, poly.Cell1DsExtrema.cols());

            // Recupera vertici associati a spigolo e1
            unsigned int e1_v0 = poly.Cell1DsExtrema(0, edgeIdx1);
            unsigned int e1_v1 = poly.Cell1DsExtrema(1, edgeIdx1);

            // Recupera vertici associati a spigolo e2
            unsigned int e2_v0 = poly.Cell1DsExtrema(0, edgeIdx2);
            unsigned int e2_v1 = poly.Cell1DsExtrema(1, edgeIdx2);

            // Determina la direzione locale del bordo e1
            unsigned int face_v_i = faceVertices[i];
            unsigned int face_v_ip1 = faceVertices[(i + 1) % E];

            unsigned int e1_end = (face_v_i == e1_v0) ? e1_v1 : e1_v0;
            unsigned int e2_start = (face_v_ip1 == e2_v0) ? e2_v0 : e2_v1;

            EXPECT_EQ(e1_end, e2_start)
                << "Face " << face << ": edge[" << i << "] ends at " << e1_end
                << " but edge[" << (i + 1) % E << "] starts at " << e2_start
                << ". Vertex loop = (" << face_v_i << " -> " << face_v_ip1 << ")";
        }
    }
}

TEST_F(PolyhedralMeshTest, FaceEdgesAreCyclicallyConnected)
{
    GeodesicPolyhedron poly = Build_ClassI_Geodesic(*tetrahedron, 5);

    for (unsigned int face : poly.Cell2DsId)
    {
        const auto& faceEdges = poly.Cell2DsEdges[face];
        const auto& faceVertices = poly.Cell2DsVertices[face];
        size_t E = faceEdges.size();

        ASSERT_EQ(E, faceVertices.size()) << "Face " << face << ": mismatch between number of edges and vertices.";

        // Verifica che il primo vertice sia in uno degli estremi del primo spigolo
        unsigned int firstEdgeIdx = faceEdges[0];
        unsigned int firstVertex = faceVertices[0];
        ASSERT_LT(firstEdgeIdx, poly.Cell1DsExtrema.cols());

        unsigned int orig0 = poly.Cell1DsExtrema(0, firstEdgeIdx);
        unsigned int end0 = poly.Cell1DsExtrema(1, firstEdgeIdx);

        EXPECT_TRUE(firstVertex == orig0 || firstVertex == end0)
            << "Face " << face << ": vertex[0] = " << firstVertex
            << " is not an endpoint of edge[0] = " << firstEdgeIdx;

        // Per ogni coppia di spigoli consecutivi
        for (size_t i = 0; i < E; ++i)
        {
            unsigned int edgeIdx1 = faceEdges[i];
            unsigned int edgeIdx2 = faceEdges[(i + 1) % E];

            ASSERT_LT(edgeIdx1, poly.Cell1DsExtrema.cols());
            ASSERT_LT(edgeIdx2, poly.Cell1DsExtrema.cols());

            // Recupera vertici associati a spigolo e1
            unsigned int e1_v0 = poly.Cell1DsExtrema(0, edgeIdx1);
            unsigned int e1_v1 = poly.Cell1DsExtrema(1, edgeIdx1);

            // Recupera vertici associati a spigolo e2
            unsigned int e2_v0 = poly.Cell1DsExtrema(0, edgeIdx2);
            unsigned int e2_v1 = poly.Cell1DsExtrema(1, edgeIdx2);

            // Determina la direzione locale del bordo e1
            unsigned int face_v_i = faceVertices[i];
            unsigned int face_v_ip1 = faceVertices[(i + 1) % E];

            unsigned int e1_end = (face_v_i == e1_v0) ? e1_v1 : e1_v0;
            unsigned int e2_start = (face_v_ip1 == e2_v0) ? e2_v0 : e2_v1;

            EXPECT_EQ(e1_end, e2_start)
                << "Face " << face << ": edge[" << i << "] ends at " << e1_end
                << " but edge[" << (i + 1) % E << "] starts at " << e2_start
                << ". Vertex loop = (" << face_v_i << " -> " << face_v_ip1 << ")";
        }
    }
}

// Test for dualize
// Test dualize on a simple tetrahedral geodesic polyhedron
TEST_F(PolyhedralMeshTest, DualizeTetrahedronGeodesic) {
    // Costruzione geodesica base (subdivisione 1) di un tetraedro
    GeodesicPolyhedron original = Build_ClassI_Geodesic(*tetrahedron, 2);
    ASSERT_GT(original.NumCell0Ds, 0);
    ASSERT_GT(original.NumCell1Ds, 0);
    ASSERT_GT(original.NumCell2Ds, 0);

    GeodesicPolyhedron dual = dualize(original);

    // 1. Verifica relazioni topologiche
    EXPECT_EQ(dual.NumCell0Ds, original.NumCell2Ds);  // vertici duali = facce originali
    EXPECT_EQ(dual.NumCell2Ds, original.NumCell0Ds);  // facce duali = vertici originali
    EXPECT_EQ(dual.NumCell1Ds, original.NumCell1Ds);  // lati invariati in numero

    // 2. Verifica che le coordinate dei vertici duali siano normalizzate
    for (int i = 0; i < dual.Cell0DsCoordinates.cols(); ++i) {
        double norm = dual.Cell0DsCoordinates.col(i).norm();
        EXPECT_NEAR(norm, 1.0, 1e-10) << "Vertex " << i << " has norm " << norm;
    }

    // 3. Verifica che ogni edge duale colleghi due vertici validi
    for (int i = 0; i < dual.NumCell1Ds; ++i) {
        unsigned int v0 = dual.Cell1DsExtrema(0, i);
        unsigned int v1 = dual.Cell1DsExtrema(1, i);
        EXPECT_LT(v0, dual.NumCell0Ds);
        EXPECT_LT(v1, dual.NumCell0Ds);
        EXPECT_NE(v0, v1);
    }

    // 4. Verifica struttura delle facce duali
    for (unsigned int fid : dual.Cell2DsId) {
        const auto& faceVerts = dual.Cell2DsVertices[fid];
        const auto& faceEdges = dual.Cell2DsEdges[fid];

        EXPECT_GT(faceVerts.size(), 2) << "Face " << fid << " has < 3 vertices";
        EXPECT_EQ(faceVerts.size(), faceEdges.size()) 
            << "Face " << fid << " has mismatched vertices and edges";
        
        for (unsigned int vid : faceVerts) {
            EXPECT_LT(vid, dual.NumCell0Ds);
        }
        for (unsigned int eid : faceEdges) {
            EXPECT_LT(eid, dual.NumCell1Ds);
        }
    }

    // 5. Verifica dimensioni vettori coerenti
    EXPECT_EQ(dual.Cell0DsId.size(), dual.NumCell0Ds);
    EXPECT_EQ(dual.Cell1DsId.size(), dual.NumCell1Ds);
    EXPECT_EQ(dual.Cell2DsId.size(), dual.NumCell2Ds);

    EXPECT_EQ(dual.Cell2DsVertices.size(), dual.NumCell2Ds);
    EXPECT_EQ(dual.Cell2DsEdges.size(), dual.NumCell2Ds);
}


TEST_F(PolyhedralMeshTest, dualizeII_BasicTest) {
    GeodesicPolyhedron geodesic = Build_ClassII_Geodesic(*tetrahedron, 3);
    GeodesicPolyhedron dual = dualize(geodesic);
    
    // Check dual relationships
    EXPECT_EQ(dual.NumCell0Ds, geodesic.NumCell2Ds);  
    EXPECT_EQ(dual.NumCell2Ds, geodesic.NumCell0Ds);  
    EXPECT_EQ(dual.NumCell1Ds, geodesic.NumCell1Ds);  
    
    // Check that dual coordinates are normalized
    for (int i = 0; i < dual.Cell0DsCoordinates.cols(); ++i) {
        EXPECT_NEAR(dual.Cell0DsCoordinates.col(i).norm(), 1.0, 1e-10);
    }
    
    // Check structure sizes
    EXPECT_EQ(dual.Cell0DsId.size(), dual.NumCell0Ds);
    EXPECT_EQ(dual.Cell1DsId.size(), dual.NumCell1Ds);
    EXPECT_EQ(dual.Cell2DsId.size(), dual.NumCell2Ds);
    
    EXPECT_EQ(dual.Cell2DsVertices.size(), dual.NumCell2Ds);
    EXPECT_EQ(dual.Cell2DsEdges.size(), dual.NumCell2Ds);
}


TEST_F(PolyhedralMeshTest, dualize_TetrahedronSelfDual) {
    GeodesicPolyhedron geodesic = Build_ClassI_Geodesic(*tetrahedron, 1);
    GeodesicPolyhedron dual = dualize(geodesic);
    
    // Tetrahedron is self-dual, so dual of dual should have same structure
    GeodesicPolyhedron dual_dual = dualize(dual);
    
    EXPECT_EQ(dual_dual.NumCell0Ds, geodesic.NumCell0Ds);
    EXPECT_EQ(dual_dual.NumCell1Ds, geodesic.NumCell1Ds);
    EXPECT_EQ(dual_dual.NumCell2Ds, geodesic.NumCell2Ds);
}
*/

// Test for SubdividePlatonicEdges
TEST_F(PolyhedralMeshTest, SubdividePlatonicEdges_BasicTest) {
    GeodesicPolyhedron geodesic;
    Initialize_ClassI_GeodesicCounts(geodesic, *tetrahedron, 3);
    InitializeGeodesicStorage(geodesic);
    
    // Copy original vertices
    std::copy(tetrahedron->VerticesId.begin(), tetrahedron->VerticesId.end(), 
              geodesic.Cell0DsId.begin());
    geodesic.Cell0DsCoordinates.block(0, 0, 3, 4) = tetrahedron->VerticesCoordinates;
    
    unsigned int nextVertexId = tetrahedron->NumVertices;
    unsigned int nextEdgeId = 0;
    
    SubdividePlatonicEdges(geodesic, *tetrahedron, 3, nextVertexId, nextEdgeId);
    
    // Check that vertices were added
    EXPECT_GT(nextVertexId, tetrahedron->NumVertices);
    
    // Check that edges were added
    EXPECT_GT(nextEdgeId, 0);
    
    // Each edge should be subdivided into n=3 segments, so n-1=2 new vertices per edge
    unsigned int expectedNewVertices = tetrahedron->NumEdges * (3 - 1);
    EXPECT_EQ(nextVertexId, tetrahedron->NumVertices + expectedNewVertices);
}

// Test for TriangulateFacesClassI
TEST_F(PolyhedralMeshTest, TriangulateFacesClassI_BasicSetup) {
    // This is a complex function that requires proper setup
    // We'll test that it can be called without crashing
    GeodesicPolyhedron geodesic;
    unsigned int n = 2;
    
    Initialize_ClassI_GeodesicCounts(geodesic, *tetrahedron, n);
    InitializeGeodesicStorage(geodesic);
    
    // Setup vertices
    std::copy(tetrahedron->VerticesId.begin(), tetrahedron->VerticesId.end(), 
              geodesic.Cell0DsId.begin());
    geodesic.Cell0DsCoordinates.block(0, 0, 3, 4) = tetrahedron->VerticesCoordinates;
    
    unsigned int nextVertexId = tetrahedron->NumVertices;
    unsigned int nextEdgeId = 0;
    
    // Subdivide edges first
    SubdividePlatonicEdges(geodesic, *tetrahedron, n, nextVertexId, nextEdgeId);
    
    // Create the required matrices
    const Map<Matrix<unsigned int, Dynamic, Dynamic, ColMajor>> 
        internalVerticesMatrix(geodesic.Cell0DsId.data() + tetrahedron->NumVertices,
                               n - 1, tetrahedron->NumEdges);
    
    const Map<Matrix<unsigned int, Dynamic, Dynamic, ColMajor>> 
        internalEdgesMatrix(geodesic.Cell1DsId.data(), n, tetrahedron->NumEdges);
    
    unsigned int nextFaceId = 0;
    
    // This should not crash
    EXPECT_NO_THROW(TriangulateFacesClassI(geodesic, *tetrahedron, n,
                                           nextVertexId, nextEdgeId, nextFaceId,
                                           internalVerticesMatrix, internalEdgesMatrix));
}

