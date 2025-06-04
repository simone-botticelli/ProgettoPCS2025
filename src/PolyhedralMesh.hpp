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

struct GeodesicCounts {
    unsigned int V;
    unsigned int E;
    unsigned int F;
};

GeodesicCounts SetGeodesicCounts_ClassI(const PlatonicSolid& solid, const unsigned int& n);

void addVertex(GeodesicPolyhedron& geodesic, unsigned int vertexId, const Vector3d& vertexCoordinates);

void addEdge(GeodesicPolyhedron& geodesic, unsigned int edgeId, unsigned int originId, unsigned int endId);

unsigned int GetorAddEdge(GeodesicPolyhedron& geodesic, unsigned int& nextEdgeId, unsigned int originId, unsigned int endId);

void addFace( 
    GeodesicPolyhedron& geodesic,
    unsigned int& nextEdgeId,
    unsigned int faceId,
    const vector<unsigned int>& facesVertices
);

// Overload per triangoli
void addFace(
    GeodesicPolyhedron& geodesic,
    unsigned int& nextEdgeId,
    unsigned int faceId,
    unsigned int v0, unsigned int v1, unsigned int v2);

GeodesicPolyhedron Build_ClassI_Geodesic(const PlatonicSolid& solid, const unsigned int n);

}