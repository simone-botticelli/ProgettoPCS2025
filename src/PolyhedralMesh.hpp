#pragma once

#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace PolyhedralLibrary {

struct GeodesicPolyhedron
{
    unsigned int NumCell0Ds;
    unsigned int NumCell1Ds;
    unsigned int NumCell2Ds;

    vector<unsigned int> Cell0DsId;
    vector<unsigned int> Cell1DsId;
    vector<unsigned int> Cell2DsId;
    unsigned int Cell3DsId;
    
    vector<bool> Cell0DsShortPath;
    vector<bool> Cell1DsShortPath;

	map<pair<unsigned int, unsigned int>, unsigned int> ExtrematoEdge;
    MatrixXd Cell0DsCoordinates;
    MatrixXi Cell1DsExtrema;
    vector<vector<unsigned int>> Cell2DsVertices;
    vector<vector<unsigned int>> Cell2DsEdges;
};

enum class PlatonicType {
    TETRAHEDRON,
    OCTAHEDRON,
    ICOSAHEDRON
};

struct PlatonicSolid 
{
	const PlatonicType Type;
	unsigned int NumVertices, NumEdges, NumFaces;
	vector<unsigned int> VerticesId, EdgesId, FacesId;
	MatrixXd VerticesCoordinates;
	MatrixXi EdgesExtrema;
	vector<vector<unsigned int>> FacesVertices;
	vector<vector<unsigned int>> FacesEdges;
	
	PlatonicSolid(PlatonicType type);
	
};

void NormalizeMatrixColumns(MatrixXd& M);

void Initialize_ClassI_GeodesicCounts(GeodesicPolyhedron& geodesic,const PlatonicSolid& solid, const unsigned int& n);

void Initialize_ClassII_GeodesicCounts(GeodesicPolyhedron& geodesic,const PlatonicSolid& solid, const unsigned int& n);

void InitializeGeodesicStorage(GeodesicPolyhedron& geodesic);

void addVertex(GeodesicPolyhedron& geodesic, unsigned int vertexId, const Vector3d& vertexCoordinates);

unsigned int GetorAddEdge(GeodesicPolyhedron& geodesic, unsigned int& nextEdgeId, unsigned int originId, unsigned int endId);

/*void addFace( 
    GeodesicPolyhedron& geodesic,
    unsigned int& nextEdgeId,
    unsigned int faceId,
    const vector<unsigned int>& facesVertices
);*/

// Overload per triangoli
void addFace(
    GeodesicPolyhedron& geodesic,
    unsigned int& nextEdgeId,
    unsigned int& nextFaceId,
    unsigned int v0, unsigned int v1, unsigned int v2);
    
void SubdividePlatonicEdges(GeodesicPolyhedron& geodesic,
                              const PlatonicSolid& solid,
                              unsigned int n,
                              unsigned int& nextVertexId, 
                              unsigned int& nextEdgeId);

VectorXd ComputePointOnTriangle(unsigned int n,
								unsigned int i,
								unsigned int j,
								const Ref<const VectorXd>& vA_coords, 
								const Ref<const VectorXd>& vB_coords,
								const Ref<const VectorXd>& vC_coords);

void TriangulateFacesClassI(GeodesicPolyhedron& geodesic,
                                const PlatonicSolid& solid,
                                unsigned int n,
                                unsigned int& nextVertexId,
                                unsigned int& nextEdgeId,
                                unsigned int& nextFaceId, 
                                const Map<Matrix<unsigned int, Dynamic, Dynamic, ColMajor>>& internalVerticesMatrix,
                                const Map<Matrix<unsigned int, Dynamic, Dynamic, ColMajor>>& internalEdgesMatrix);

GeodesicPolyhedron Build_ClassI_Geodesic(const PlatonicSolid& solid, const unsigned int n);

void ComputeShortestPath(GeodesicPolyhedron& mesh, unsigned int source, unsigned int target);

GeodesicPolyhedron dualize(const GeodesicPolyhedron& poly);

}