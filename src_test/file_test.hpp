#include <gtest/gtest.h>
#include "PolyhedralMesh.hpp"
#include <fstream>
#include <cmath>

using namespace PolyhedralLibrary;

class PolyhedralMeshTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Setup comune per i test
    }
    
    void TearDown() override {
        // Cleanup dopo i test
    }
};

// Test 1: Verifica costruzione Tetraedro
TEST_F(PolyhedralMeshTest, TetrahedronConstruction) {
    PlatonicSolid tetra(PlatonicType::TETRAHEDRON);
    
    EXPECT_EQ(tetra.NumVertices, 4);
    EXPECT_EQ(tetra.NumEdges, 6);
    EXPECT_EQ(tetra.NumFaces, 4);
    
    // Verifica che gli ID siano sequenziali da 0
    for(unsigned int i = 0; i < tetra.NumVertices; i++) {
        EXPECT_EQ(tetra.VerticesId[i], i);
    }
    
    // Verifica dimensioni coordinate
    EXPECT_EQ(tetra.VerticesCoordinates.rows(), 3);
    EXPECT_EQ(tetra.VerticesCoordinates.cols(), 4);
}

// Test 2: Verifica costruzione Ottaedro
TEST_F(PolyhedralMeshTest, OctahedronConstruction) {
    PlatonicSolid octa(PlatonicType::OCTAHEDRON);
    
    EXPECT_EQ(octa.NumVertices, 6);
    EXPECT_EQ(octa.NumEdges, 12);
    EXPECT_EQ(octa.NumFaces, 8);
    
    // Controllo che le facce siano triangolari
    for(unsigned int f = 0; f < octa.NumFaces; f++) {
        EXPECT_EQ(octa.FacesVertices[f].size(), 3);
        EXPECT_EQ(octa.FacesEdges[f].size(), 3);
    }
}

// Test 3: Verifica costruzione Icosaedro
TEST_F(PolyhedralMeshTest, IcosahedronConstruction) {
    PlatonicSolid icosa(PlatonicType::ICOSAHEDRON);
    
    EXPECT_EQ(icosa.NumVertices, 12);
    EXPECT_EQ(icosa.NumEdges, 30);
    EXPECT_EQ(icosa.NumFaces, 20);
    
    // Verifica che tutti i vertici abbiano coordinate valide
    for(int i = 0; i < icosa.VerticesCoordinates.cols(); i++) {
        Vector3d vertex = icosa.VerticesCoordinates.col(i);
        EXPECT_FALSE(std::isnan(vertex.norm()));
        EXPECT_GT(vertex.norm(), 0.0);
    }
}

// Test 4: Formula di Eulero
TEST_F(PolyhedralMeshTest, EulerFormula) {
    vector<PlatonicType> types = {PlatonicType::TETRAHEDRON, 
                                  PlatonicType::OCTAHEDRON, 
                                  PlatonicType::ICOSAHEDRON};
    
    for(auto type : types) {
        PlatonicSolid solid(type);
        // V - E + F = 2 per poliedri convessi
        int euler = solid.NumVertices - solid.NumEdges + solid.NumFaces;
        EXPECT_EQ(euler, 2) << "Formula di Eulero non rispettata per tipo: " << (int)type;
    }
}

// Test 5: Conteggi geodetici Classe I
TEST_F(PolyhedralMeshTest, GeodesicCountsClassI) {
    // Test per n = 1
    GeodesicPolyhedron g1;
    PlatonicSolid tetra(PlatonicType::TETRAHEDRON);
    Initialize_ClassI_GeodesicCounts(g1, tetra, 1);
    EXPECT_EQ(g1.NumCell0Ds, 4);   // 2*1^2 + 2 = 4
    EXPECT_EQ(g1.NumCell1Ds, 6);   // 6*1^2 = 6
    EXPECT_EQ(g1.NumCell2Ds, 4);   // 4*1^2 = 4

    // Test per n = 2
    GeodesicPolyhedron g2;
    Initialize_ClassI_GeodesicCounts(g2, tetra, 2);
    EXPECT_EQ(g2.NumCell0Ds, 10);  // 2*4 + 2 = 10
    EXPECT_EQ(g2.NumCell1Ds, 24);  // 6*4 = 24
    EXPECT_EQ(g2.NumCell2Ds, 16);  // 4*4 = 16
}


// Test 6: Costruzione geodetica base
TEST_F(PolyhedralMeshTest, GeodesicBasicConstruction) {
    PlatonicSolid tetra(PlatonicType::TETRAHEDRON);
    GeodesicPolyhedron geo = Build_ClassI_Geodesic(tetra, 1);
    
    // Per n=1 dovrebbe coincidere con il tetraedro originale
    EXPECT_EQ(geo.NumCell0Ds, tetra.NumVertices);
    EXPECT_EQ(geo.Cell0DsId.size(), tetra.NumVertices);
    
    // Verifica che le coordinate siano copiate correttamente
    for(unsigned int i = 0; i < tetra.NumVertices; i++) {
        VectorXd orig = tetra.VerticesCoordinates.col(i);
        VectorXd geo_v = geo.Cell0DsCoordinates.col(i);
        EXPECT_NEAR((orig - geo_v).norm(), 0.0, 1e-10);
    }
}

// Test 7: Verifica input errati
TEST_F(PolyhedralMeshTest, InvalidInputHandling) {
    PlatonicSolid tetra(PlatonicType::TETRAHEDRON);
    GeodesicPolyhedron geo = Build_ClassI_Geodesic(tetra, 2);
    
    // Test con indici vertice non validi per shortest path
    testing::internal::CaptureStderr();
    ComputeShortestPath(geo, 999, 0);  // indice non valido
    string error_output = testing::internal::GetCapturedStderr();
    EXPECT_FALSE(error_output.empty());
}

// Test 8: Proprietà geometriche - vertici sulla sfera unitaria
TEST_F(PolyhedralMeshTest, VerticesOnUnitSphere) {
    PlatonicSolid icosa(PlatonicType::ICOSAHEDRON);
    
    // Verifica che i vertici siano normalizzati (approssimativamente sulla sfera unitaria)
    for(int i = 0; i < icosa.VerticesCoordinates.cols(); i++) {
        Vector3d vertex = icosa.VerticesCoordinates.col(i);
        double norm = vertex.norm();
        // Tolleranza per errori numerici
        EXPECT_NEAR(norm, 1.0, 1e-1) << "Vertice " << i << " non sulla sfera unitaria";
    }
}

// Test 9: Struttura dati corretta per celle 1D
TEST_F(PolyhedralMeshTest, EdgeStructureValid) {
    PlatonicSolid tetra(PlatonicType::TETRAHEDRON);
    
    // Verifica che ogni edge abbia esattamente 2 estremi
    EXPECT_EQ(tetra.EdgesExtrema.rows(), 2);
    EXPECT_EQ(tetra.EdgesExtrema.cols(), tetra.NumEdges);
    
    // Verifica che gli estremi siano indici validi di vertici
    for(int e = 0; e < tetra.NumEdges; e++) {
        int v1 = tetra.EdgesExtrema(0, e);
        int v2 = tetra.EdgesExtrema(1, e);
        
        EXPECT_GE(v1, 0);
        EXPECT_LT(v1, (int)tetra.NumVertices);
        EXPECT_GE(v2, 0);
        EXPECT_LT(v2, (int)tetra.NumVertices);
        EXPECT_NE(v1, v2); // Un edge non può collegare un vertice a se stesso
    }
}

// Test 10: Consistenza delle facce
TEST_F(PolyhedralMeshTest, FaceConsistency) {
    PlatonicSolid octa(PlatonicType::OCTAHEDRON);
    
    for(unsigned int f = 0; f < octa.NumFaces; f++) {
        vector<unsigned int>& vertices = octa.FacesVertices[f];
        vector<unsigned int>& edges = octa.FacesEdges[f];
        
        // Numero di vertici = numero di edge per facce triangolari
        EXPECT_EQ(vertices.size(), edges.size());
        
        // Verifica che tutti gli indici siano validi
        for(unsigned int v : vertices) {
            EXPECT_LT(v, octa.NumVertices);
        }
        for(unsigned int e : edges) {
            EXPECT_LT(e, octa.NumEdges);
        }
    }
}

// Test 12: Shortest path base
TEST_F(PolyhedralMeshTest, ShortestPathBasic) {
    PlatonicSolid tetra(PlatonicType::TETRAHEDRON);
    GeodesicPolyhedron geo = Build_ClassI_Geodesic(tetra, 2);
    
    // Inizializza i flag di shortest path
    geo.Cell0DsShortPath.resize(geo.NumCell0Ds, false);
    geo.Cell1DsShortPath.resize(geo.NumCell1Ds, false);
    
    // Test con percorso valido
    testing::internal::CaptureStdout();
    ComputeShortestPath(geo, 0, 1);
    string output = testing::internal::GetCapturedStdout();
    
    // Verifica che sia stato trovato un percorso
    EXPECT_TRUE(output.find("Cammino minimo trovato") != string::npos);
}

// Test 13: Dualization base  
TEST_F(PolyhedralMeshTest, DualizationBasic) {
    PlatonicSolid tetra(PlatonicType::TETRAHEDRON);
    GeodesicPolyhedron geo = Build_ClassI_Geodesic(tetra, 2);
    
    // Assicura che la struttura sia completa per la dualizzazione
    geo.Cell2DsVertices = tetra.FacesVertices;
    geo.Cell2DsEdges = tetra.FacesEdges;
    
    GeodesicPolyhedron dual = dualize(geo);
    
    // Verifica proprietà duali
    EXPECT_EQ(dual.NumCell0Ds, geo.NumCell2Ds); // vertici dual = facce originali
    EXPECT_EQ(dual.NumCell2Ds, geo.NumCell0Ds); // facce dual = vertici originali
    EXPECT_EQ(dual.NumCell1Ds, geo.NumCell1Ds); // edges rimangono uguali
}
    
