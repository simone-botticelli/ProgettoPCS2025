#include <iostream>
#include <cassert>
#include <queue>
#include <limits>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <span>

#include "PolyhedralMesh.hpp"
#include "Eigen/Eigen"


using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary {

PlatonicSolid::PlatonicSolid(PlatonicType type) : Type(type)
{	
	switch(type)
	{
		case PlatonicType::TETRAHEDRON:{
			NumVertices = 4;
			NumEdges = 6;
			NumFaces = 4;
			
			VerticesId = {0, 1, 2, 3};
			EdgesId    = {0, 1, 2, 3, 4, 5};
			FacesId    = {0, 1, 2, 3};
			
			VerticesCoordinates.resize(3, NumVertices);
			VerticesCoordinates << 
				 1, -1, -1,  1,
				 1, -1,  1, -1,
				 1,  1, -1, -1;
			
			EdgesExtrema.resize(2, NumEdges);
			EdgesExtrema << 
				0, 0, 0, 1, 1, 2,
				1, 2, 3, 2, 3, 3;
			
			FacesVertices = {
				{0, 1, 2},   
				{0, 1, 3},   
				{0, 2, 3},   
				{1, 2, 3}    
			};
	
			FacesEdges = {
				{0, 3, 1},
				{0, 4, 2},  
				{1, 5, 2},  
				{3, 5, 4}    
			};
			break;
		}		
		case PlatonicType::OCTAHEDRON:{
			NumVertices = 6;
			NumEdges = 12;
			NumFaces = 8;
	
			VerticesId = {0, 1, 2, 3, 4, 5};
			EdgesId = {0, 1, 2, 3, 4, 5, 6, 7 ,8, 9 ,10, 11};
			FacesId = {0, 1, 2, 3, 4, 5, 6, 7};
	
			VerticesCoordinates.resize(3, NumVertices);
			VerticesCoordinates <<
				0,  1,  0, -1,  0,  0,
				0,  0,  1,  0, -1,  0,
				1,  0,  0,  0,  0, -1;
	
			EdgesExtrema.resize(2, NumEdges);
			EdgesExtrema <<
			0, 0, 0, 0, 1, 2, 3, 1, 1, 2, 3, 4,
			1, 2, 3, 4, 2, 3, 4, 4, 5, 5, 5, 5;
	
			FacesVertices = {
				{0, 1, 2},
				{0, 2, 3},
				{0, 3, 4},
				{0, 1, 4},
				{1, 2, 5},  
				{2, 3, 5},
				{3, 4, 5},
				{1, 4, 5}
			};
	
			FacesEdges = {
				{0, 4, 1},
				{1, 5, 2},
				{2, 6, 3},
				{0, 7, 3},
				{4, 9, 8},
				{5, 10, 9},
				{6, 11, 10},
				{7, 11, 8}
			};
			break;
		}		
		case PlatonicType::ICOSAHEDRON:{
			const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
			
			NumVertices = 12;
			NumEdges = 30;
			NumFaces = 20;
			
			VerticesId = {0, 1, 2, 3, 4, 5, 6, 7 ,8, 9 ,10, 11};
			EdgesId = {0, 1, 2, 3, 4, 5, 6, 7 ,8, 9 ,10, 11,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19, 20, 21, 22, 23, 24 ,25, 26, 27, 28, 29};
			FacesId = {0, 1, 2, 3, 4, 5, 6, 7 ,8, 9 ,10, 11,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19};
	
			VerticesCoordinates.resize(3, NumVertices);
			VerticesCoordinates <<
				-1,  1, -1,  1,  0,  0,  0,  0,  phi,  phi, -phi, -phi,
				 phi, phi, -phi, -phi, -1,  1, -1,  1,  0,   0,   0,   0,
				 0,   0,   0,   0,  phi, phi, -phi, -phi, -1,   1,  -1,   1;
	
			EdgesExtrema.resize(2, NumEdges);
			EdgesExtrema <<
				 0, 5, 0, 1, 0, 1, 0, 7, 0,10, 5, 1, 4, 4, 2, 2, 6, 6, 1, 7, 3, 4, 3, 2, 2, 2, 3, 6, 3, 8, 
				11,11, 5, 5, 1, 7, 7,10,10,11, 9, 9,11, 5,10,11, 7,10, 8, 8, 9, 9, 4, 4, 3, 6, 6, 8, 8, 9;
	
			FacesVertices = {
				{0, 5, 11}, {0, 1, 5}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
				{1, 5, 9}, {4, 5, 11}, {2, 10, 11}, {6, 7, 10}, {1, 7, 8},
				{3, 4, 9}, {2, 3, 4}, {2, 3, 6}, {3, 6, 8}, {3, 8, 9},
				{4, 5, 9}, {2, 4, 11}, {2, 6, 10}, {6, 7, 8}, {1, 8, 9}};
			
			FacesEdges = {
				{2, 1, 0},     
				{4, 3, 2},    
				{4, 5, 6}, 
				{6, 7, 8},
				{8, 9, 0},  
				{3, 10, 11},  
				{13, 1, 12},  
				{14, 9, 15},  
				{16, 7, 17},   
				{5, 19, 18}, 
				{22, 21, 20},   
				{24, 22, 23},   
				{24, 26, 25},   
				{26, 27, 28},  
				{28, 29, 20},  
				{13, 10, 21},   
				{23, 12, 15},   
				{25, 17, 14},   
				{16, 19, 27},  
				{18, 29, 11}    
			};
			break;
		}
		default:
			cerr << "PlatonicSolid constructor is defined only for Tetrahedron, Octahedron, Icosahedron" << endl;
			break;
		}
		VerticesEdges.resize(NumVertices);
			
			// Per ogni spigolo, aggiungilo ai due vertici estremi
			for (unsigned int edgeId = 0; edgeId < NumEdges; ++edgeId) {
				unsigned int v0 = EdgesExtrema(0, edgeId);
				unsigned int v1 = EdgesExtrema(1, edgeId);
				VerticesEdges[v0].push_back(edgeId);
				VerticesEdges[v1].push_back(edgeId);
				}
};

void Initialize_ClassI_GeodesicCounts(GeodesicPolyhedron& geodesic,const PlatonicSolid& solid, const unsigned int& n)
{
	unsigned int T = n*n;
	
	switch(solid.Type) 
	{
		case PlatonicType::TETRAHEDRON:
			geodesic.NumCell0Ds = 2*T + 2;
			geodesic.NumCell1Ds = 6*T;
			geodesic.NumCell2Ds = 4*T;
				break;
			
		case PlatonicType::OCTAHEDRON:
			geodesic.NumCell0Ds = 4*T + 2;
			geodesic.NumCell1Ds = 12*T;
			geodesic.NumCell2Ds = 8*T;
				break;
		
		case PlatonicType::ICOSAHEDRON:
			geodesic.NumCell0Ds = 10*T + 2;
			geodesic.NumCell1Ds = 30*T;
			geodesic.NumCell2Ds = 20*T;
				break;
			
		default:
			cerr << "PlatonicType not supported" << endl;
			throw std::runtime_error("Unsupported PlatonicType in Initialize_ClassI_GeodesicCounts");
	}
}

void Initialize_ClassII_GeodesicCounts(GeodesicPolyhedron& geodesic,const PlatonicSolid& solid, const unsigned int& n)
{
	unsigned int T = n*n;
	
	switch(solid.Type) 
	{
		case PlatonicType::TETRAHEDRON:
			geodesic.NumCell0Ds = 6*(T + n) + 2;
			geodesic.NumCell1Ds = 18*(T + n);
			geodesic.NumCell2Ds = 12*(T + n);
				break;
			
		case PlatonicType::OCTAHEDRON:
			geodesic.NumCell0Ds = 12*(T + n) + 2;
			geodesic.NumCell1Ds = 36*(T + n);
			geodesic.NumCell2Ds = 24*(T + n);
				break;
		
		case PlatonicType::ICOSAHEDRON:
			geodesic.NumCell0Ds = 30*(T + n) + 2;
			geodesic.NumCell1Ds = 90*(T + n);
			geodesic.NumCell2Ds = 60*(T + n);
				break;
			
		default:
			cerr << "PlatonicType not supported" << endl;
			throw std::runtime_error("Unsupported PlatonicType in Initialize_ClassII_GeodesicCounts");
	}
}

void InitializeGeodesicStorage(GeodesicPolyhedron& geodesic) {
    geodesic.Cell0DsId.resize(geodesic.NumCell0Ds);
    geodesic.Cell1DsId.resize(geodesic.NumCell1Ds);
    geodesic.Cell2DsId.resize(geodesic.NumCell2Ds);
    geodesic.Cell0DsShortPath.resize(geodesic.NumCell0Ds, false);
    geodesic.Cell1DsShortPath.resize(geodesic.NumCell1Ds, false);

    geodesic.Cell0DsCoordinates.resize(3, geodesic.NumCell0Ds);
    geodesic.Cell1DsExtrema.resize(2, geodesic.NumCell1Ds);
    
    geodesic.VertextoEdges.resize(geodesic.NumCell0Ds);

    geodesic.Cell2DsVertices.resize(geodesic.NumCell2Ds);
    geodesic.Cell2DsEdges.resize(geodesic.NumCell2Ds);
    
    geodesic.Cell3DsId = 0;
}



void addVertex(GeodesicPolyhedron& geodesic, unsigned int vertexId, const VectorXd& vertexCoordinates)
{
	assert(vertexId < geodesic.Cell0DsCoordinates.cols() && "vertexId fuori range");
	geodesic.Cell0DsId[vertexId] = vertexId;
	geodesic.Cell0DsCoordinates.col(vertexId) = vertexCoordinates;
}



unsigned int GetorAddEdge(GeodesicPolyhedron& geodesic, unsigned int& nextEdgeId,const unsigned int originId,const unsigned int endId)
{
	
        // Ordiniamo i vertici per evitare duplicati tipo (3,5) vs (5,3)
        auto& edgeMap = geodesic.ExtrematoEdge;
        auto key = std::minmax(originId, endId); // restituisce una pair già ordinata
        auto it = edgeMap.find(key);
        if (it != edgeMap.end()) {
        	return it->second; // già presente
        } else {
            int newId = nextEdgeId++;
            edgeMap[key] = newId;
            geodesic.Cell1DsId[newId] = newId;
            geodesic.Cell1DsExtrema(0, newId) = originId;
			geodesic.Cell1DsExtrema(1, newId) = endId;
			geodesic.VertextoEdges[originId].push_back(newId);
			geodesic.VertextoEdges[endId].push_back(newId);
			return newId;
            }
       
}

void addEdge(GeodesicPolyhedron& geodesic, unsigned int newId, unsigned int originId, unsigned int endId)
{
	auto& edgeMap = geodesic.ExtrematoEdge;
	auto key = std::minmax(originId, endId);
	edgeMap[key] = newId;
	geodesic.Cell1DsId[newId] = newId;
	geodesic.Cell1DsExtrema(0, newId) = originId;
	geodesic.Cell1DsExtrema(1, newId) = endId;
	geodesic.VertextoEdges[originId].push_back(newId);
	geodesic.VertextoEdges[endId].push_back(newId);
}

void addFace(
    GeodesicPolyhedron& geodesic,
    unsigned int& nextEdgeId,
    unsigned int& nextFaceId,
    unsigned int v0, unsigned int v1, unsigned int v2)
{
    unsigned int newId = nextFaceId;
    
    vector<unsigned int> facesEdges(3);
    unsigned int edge1_Id = GetorAddEdge(geodesic, nextEdgeId, v0, v1);
    facesEdges[0] = edge1_Id;
    unsigned int edge2_Id = GetorAddEdge(geodesic, nextEdgeId, v1, v2);
    facesEdges[1] = edge2_Id;
    unsigned int edge3_Id = GetorAddEdge(geodesic, nextEdgeId, v2, v0);
    facesEdges[2] = edge3_Id;
    
    
    geodesic.Cell2DsVertices[newId] = {v0, v1, v2};
    facesEdges = {edge1_Id, edge2_Id, edge3_Id};
    geodesic.Cell2DsEdges[newId] = facesEdges;
    geodesic.Cell2DsId[newId] = newId;
    nextFaceId++;
}


void SubdividePlatonicEdges(GeodesicPolyhedron& geodesic,
                              const PlatonicSolid& solid,
                              unsigned int n,
                              unsigned int& nextVertexId, 
                              unsigned int& nextEdgeId)
{
	unsigned int oldEdgeId = nextEdgeId;
		
		for(int e = 0; e < solid.NumEdges; e++){

			int v1 = solid.EdgesExtrema(0, e);
			int v2   = solid.EdgesExtrema(1, e);
			
			auto [v_start, v_end] = minmax(v1, v2);

			// ottengo le coordinate dei vertici estremi
			VectorXd start = solid.VerticesCoordinates.col(v_start);
			VectorXd end   = solid.VerticesCoordinates.col(v_end);
			
			// Gestisco separatamente il primo lato piccolo
				
				// creo le coordinate del primo vertice
				double t = double(1) / n;
				VectorXd new_point = (1-t)*start + t*end;
				
				// salvo il primo vertice
				addVertex(geodesic, nextVertexId, new_point);
				
				// salvo il primo lato piccolo 
				GetorAddEdge(geodesic, nextEdgeId, v_start, nextVertexId);
				
				unsigned int oldVertexId = nextVertexId++; // assegna e poi aumenta

			// Creo gli altri vertici e lati sul lato principale
			
				for (int k = 2; k < n; ++k)  // k va da 2 a n-1
				{
					
					// creo le coordinate
					double t = double(k) / double(n);
					new_point = (1-t)*start + t*end;
					
					// salvo il vertice
					addVertex(geodesic, nextVertexId, new_point);
					
					// salvo gli estremi dei lati in Cell1DsExtrema
					GetorAddEdge(geodesic, nextEdgeId, oldVertexId, nextVertexId);
					
					oldVertexId = nextVertexId++;
					
				}
			
			// Collego l'ultimo lato piccolo
				GetorAddEdge(geodesic, nextEdgeId, oldVertexId, v_end);
			}
}


VectorXd ComputePointOnTriangle(unsigned int n,
								unsigned int i,
								unsigned int j,
								const Ref<const VectorXd>& vA_coords, 
								const Ref<const VectorXd>& vB_coords,
								const Ref<const VectorXd>& vC_coords)
{
	double alpha = double(n - i - j) / n;
	double beta  = double(i) / n;
	double gamma = double(j) / n;
	
	return alpha * vA_coords + beta * vB_coords + gamma * vC_coords;
}

void TriangulateFacesClassI(GeodesicPolyhedron& geodesic,
                                const PlatonicSolid& solid,
                                unsigned int n,
                                unsigned int& nextVertexId,
                                unsigned int& nextEdgeId,
                                unsigned int& nextFaceId, 
                                const Map<Matrix<unsigned int, Dynamic, Dynamic, ColMajor>>& internalVerticesMatrix,
                                const Map<Matrix<unsigned int, Dynamic, Dynamic, ColMajor>>& internalEdgesMatrix)
{
		// per ogni faccia del solido platonico:
		for(int f = 0; f < solid.NumFaces; f++)
		{
			// creo gli alias di id_vertici principali (A, B, C) ! CONTROLLA CHE A,B,C SIANO SALVATI IN ORDINE CRESCENTE
			const unsigned int& A = solid.FacesVertices[f][0];
			const unsigned int& B = solid.FacesVertices[f][1];
			const unsigned int& C = solid.FacesVertices[f][2];
			bool inOrder = (A < B) && (B < C); // questo controllo lo possiamo poi mettere in un test 
			if (!inOrder) {
				std::cerr << "Errore: i vertici della faccia " << f << " non sono in ordine crescente" << std::endl;
			}
			
			// creo alias per le coordinate dei vertici principali
			const Ref<const VectorXd> vA = solid.VerticesCoordinates.col(A);
			const Ref<const VectorXd> vB = solid.VerticesCoordinates.col(B);
			const Ref<const VectorXd> vC = solid.VerticesCoordinates.col(C);
			
			// creo gli alias dei vettori di id_vertici interni ai lati principali (AB, BC, AC)
			const unsigned int& edge_0 = solid.FacesEdges[f][0];
			const unsigned int& edge_1 = solid.FacesEdges[f][1];
			const unsigned int& edge_2 = solid.FacesEdges[f][2];
			
			const auto& AB = internalVerticesMatrix.col(edge_0);
			const auto& BC = internalVerticesMatrix.col(edge_1);
			const auto& AC = internalVerticesMatrix.col(edge_2);
			
			// creo gli alias al vettore dei lati interni ai lati principali (AB_edges, BC_edges, AC_edges)
			const auto& AB_edges = internalEdgesMatrix.col(edge_0);
			const auto& BC_edges = internalEdgesMatrix.col(edge_1);
			const auto& AC_edges = internalEdgesMatrix.col(edge_2);
					
			// caso n = 2
			if(n==2){
				addFace(geodesic, nextEdgeId, nextFaceId, 
							A, AB[0], AC[0]);
				addFace(geodesic, nextEdgeId, nextFaceId, 
							B, AB[0], BC[0]);
				addFace(geodesic, nextEdgeId, nextFaceId, 
							C, AC[0], BC[0]);
				addFace(geodesic, nextEdgeId, nextFaceId, 
							BC[0], AB[0], AC[0]);
							
				continue;
				}
							
			// caso n > 2 
				
				// passo 1: salvo la prima faccia 
					addFace(geodesic, nextEdgeId, nextFaceId, 
							A, AB[0], AC[0]);
				
				// passo 2: creo il primo vertice interno
					int i = 1;
					int j = 1;
					VectorXd new_point = ComputePointOnTriangle(n, i, j, vA, vB, vC);
					addVertex(geodesic, nextVertexId, new_point);
					
					// salvo le tre facce tra il primo vertice interno ed A
						addFace(geodesic, nextEdgeId, nextFaceId,
								AB[0], AB[1], nextVertexId);
						addFace(geodesic, nextEdgeId, nextFaceId,
								nextVertexId, AC[0], AB[0]);
						addFace(geodesic, nextEdgeId, nextFaceId,
								nextVertexId, AC[1], AC[0]);
					nextVertexId++;
										 
				// passo 3: completo la triangolazione per righe (tranne l'ultima)
				for (int s = 3; s < n; s++)
				{
					// Creo il primo vertice interno (j = 1)
					j = 1;
					i = s - j;
					new_point = ComputePointOnTriangle(n, i, j, vA, vB, vC);
					addVertex(geodesic, nextVertexId, new_point);
						
						// salvo le due facce adiacenti ad AB
						addFace(geodesic, nextEdgeId, nextFaceId,
								AB[i - 1], AB[i], nextVertexId);
						addFace(geodesic, nextEdgeId, nextFaceId,
								AB[i - 1], nextVertexId, nextVertexId - i + 1);

					nextVertexId++;
						
					// Creo il secondo vertice interno
					j = 2;
					i = s - j;
					new_point = ComputePointOnTriangle(n, i, j, vA, vB, vC);
					addVertex(geodesic, nextVertexId, new_point); 
						
						
						// Creo la successiva lungo la riga 
						addFace(geodesic, nextEdgeId, nextFaceId,
								nextVertexId - 1, nextVertexId, nextVertexId - s + 1);
						
					nextVertexId++;
					
					// Per ogni colonna (da 2 alla penultima)
					for (j = 3; j < s; j++)
					{
						
						i = s - j;

						// creo la faccia alla colonna j 
						addFace(geodesic, nextEdgeId, nextFaceId,
								nextVertexId - 1, nextVertexId - s + 1, nextVertexId - s);
						

						// aggiungo il vertice della colonna j + 1 
						new_point = ComputePointOnTriangle(n, i, j, vA, vB, vC);
						addVertex(geodesic, nextVertexId, new_point);
						
						// creo la faccia tra la colonna j e j + 1
						addFace(geodesic, nextEdgeId, nextFaceId,
								nextVertexId - s + 1, nextVertexId - 1, nextVertexId);
						
						nextVertexId++;
					}
						
					// ultima colonna
						
						// collego l'ultimo vertice ad AC[j - 1] e creo la faccia
						
						addFace(geodesic, nextEdgeId, nextFaceId,
								AC[s-2], nextVertexId - s , nextVertexId - 1);
						
						// collego l'ultimo vertice ad AC[j] e creo la faccia
						
						addFace(geodesic, nextEdgeId, nextFaceId,
								AC[s-1], AC[s-2], nextVertexId - 1);
						
				}
				// passo 4: completo l'ultima riga
					 
					// collego le prime tre facce 
					addFace(geodesic, nextEdgeId, nextFaceId,
								AB[n-2], B, BC[0]);
					addFace(geodesic, nextEdgeId, nextFaceId,
								AB[n-2], BC[0], nextVertexId - n + 2);
					addFace(geodesic, nextEdgeId, nextFaceId,
								nextVertexId - n + 2, BC[0], BC[1]);
								
					// per ogni vertice sul lato BC (dal secondo al penultimo[1:n-3]) 
					for(int k = 1; k < n - 2; k++)
					{
						// salvo le facce interne all'ultima riga
						int vertex_under = nextVertexId - n + 2 + k;
						addFace(geodesic, nextEdgeId, nextFaceId,
								vertex_under - 1, BC[k], vertex_under);
						addFace(geodesic, nextEdgeId, nextFaceId,
								BC[k], BC[k + 1], vertex_under);
						
					}
						
					// salvo le ultime 2 facce
						
						addFace(geodesic, nextEdgeId, nextFaceId,
								BC[n-2], AC[n-2], nextVertexId - 1);
						addFace(geodesic, nextEdgeId, nextFaceId,
								BC[n-2], C, AC[n-2]);
	
		}  // fine ciclo sulle facce 
}

void NormalizeMatrixColumns(MatrixXd& M)
{
		for (int j = 0; j < M.cols(); ++j) {
			double norm = M.col(j).norm();
			if (norm > 1e-12) {
				M.col(j) /= norm;
			}
		}

}

GeodesicPolyhedron Build_ClassI_Solid(const PlatonicSolid& solid, const unsigned int n)
{
	GeodesicPolyhedron geodesic;

	// 1. Inizializzo NumCell0Ds, NumCell1Ds, NumCell2Ds
	
		Initialize_ClassI_GeodesicCounts(geodesic, solid, n);
	
	// 2. Inizializzo le dimensioni delle strutture dati
		
		InitializeGeodesicStorage(geodesic);

	// 3. Salvo i dati dei vertici del solido platonico (id e coordinate) in  Cell0DsCoordinates e in Cell0DsId
		
		copy(solid.VerticesId.begin(), solid.VerticesId.end(), geodesic.Cell0DsId.begin());
		geodesic.Cell0DsCoordinates.block(0, 0, solid.VerticesCoordinates.rows(), solid.VerticesCoordinates.cols()) = solid.VerticesCoordinates;
	
		// se n = 1 copio struttura solido platonico
		
		if (n == 1)
		{
			geodesic.Cell1DsId = solid.EdgesId;
			geodesic.Cell2DsId = solid.FacesId;
			geodesic.Cell1DsExtrema = solid.EdgesExtrema;
			geodesic.Cell2DsVertices = solid.FacesVertices;
			geodesic.Cell2DsEdges = solid.FacesEdges;
			geodesic.VertextoEdges = solid.VerticesEdges;
			geodesic.Cell3DsId = 0;
			
			NormalizeMatrixColumns(geodesic.Cell0DsCoordinates);
			
			return geodesic;
		}
		unsigned int nextVertexId = solid.NumVertices;
		unsigned int nextEdgeId = 0;
	// 4.  SUDDIVISIONE LATI PRINCIPALI:
		SubdividePlatonicEdges(geodesic,
                              solid,
                              n,
                              nextVertexId, 
                              nextEdgeId);
		
		
			//assert(geodesic.Cell0DsId.size() >= solid.NumVertices + solid.NumEdges * (n - 1));
			
		// Creo una matrice che nell'i-esima colonna ha gli id dei vertici interni al lato i-esimo del solido
			const Map<Matrix<unsigned int, Dynamic, Dynamic, ColMajor>> 
				internalVerticesMatrix(geodesic.Cell0DsId.data() + solid.NumVertices, // punto di partenza nel vector
									   n - 1, // elementi per ogni colonna 
									   solid.NumEdges); // num colonne
	
		// Creo una matrice che nell'i-esima colonna ha gli id dei lati interni al lato i-esimo del solido
			const Map<Matrix<unsigned int, Dynamic, Dynamic, ColMajor>> 
				internalEdgesMatrix(geodesic.Cell1DsId.data(), 
									n, 
									solid.NumEdges);
			
	// 6. TRIANGOLAZIONE
		unsigned int nextFaceId = 0;
		TriangulateFacesClassI(geodesic, solid, n,
							   nextVertexId, 
							   nextEdgeId, 
							   nextFaceId,
							   internalVerticesMatrix,
							   internalEdgesMatrix);
	return geodesic;
	
} // fine Build_ClassI_Solid


GeodesicPolyhedron Build_ClassI_Geodesic(const PlatonicSolid& solid, const unsigned int n)
{
	GeodesicPolyhedron geodesic = Build_ClassI_Solid(solid, n);
	NormalizeMatrixColumns(geodesic.Cell0DsCoordinates);
	return geodesic;
}

GeodesicPolyhedron Build_ClassII_Geodesic(const PlatonicSolid& solid, const unsigned int n)
{
	GeodesicPolyhedron GeoClassI = Build_ClassI_Solid(solid, n);
	GeodesicPolyhedron GeoClassII;
	// inizializzazione 
	Initialize_ClassII_GeodesicCounts(GeoClassII, solid, n);
	InitializeGeodesicStorage(GeoClassII);
	// copio i vertici
	copy(GeoClassI.Cell0DsId.begin(), GeoClassI.Cell0DsId.end(), GeoClassII.Cell0DsId.begin());
	GeoClassII.Cell0DsCoordinates.block(0, 0, GeoClassI.Cell0DsCoordinates.rows(), GeoClassI.Cell0DsCoordinates.cols()) = GeoClassI.Cell0DsCoordinates;
	unsigned int NumPrincipalEdges = n * solid.NumEdges;
	unsigned int nextVertexId = GeoClassI.NumCell0Ds;
	
	// per ogni lato principale mi calcolo il punto medio
	int startId;
	int endId;
	VectorXd new_point;
	VectorXd start;
	VectorXd end;
	for(unsigned int e = 0; e < NumPrincipalEdges; ++e)
	{
		
		startId = GeoClassI.Cell1DsExtrema(0,e);
		endId = GeoClassI.Cell1DsExtrema(1,e);
		start = GeoClassI.Cell0DsCoordinates.col(startId);
		end   = GeoClassI.Cell0DsCoordinates.col(endId);
		new_point = 0.5 * start + 0.5 * end;
		
		addVertex(GeoClassII, nextVertexId++, new_point);
		
		
	}
	span<unsigned int> EdgesMidpoints(GeoClassII.Cell0DsId.data() + GeoClassI.NumCell0Ds, NumPrincipalEdges);

	
	// per ogni faccia di GeoClassI mi calcolo il punto medio e costruisco la matrice di aciacenza edgetoFaces

	vector<vector<unsigned int>> edgeToFaces;
	edgeToFaces.resize(GeoClassI.NumCell1Ds);
	
	int A;
	int B;
	int C;
	
	VectorXd vA;
	VectorXd vB;
	VectorXd vC;
	VectorXd midpoint;
	
	for(unsigned int face : GeoClassI.Cell2DsId)
	{
		// punto medio
		A = GeoClassI.Cell2DsVertices[face][0];
		B = GeoClassI.Cell2DsVertices[face][1];
		C = GeoClassI.Cell2DsVertices[face][2];
		
		vA = GeoClassI.Cell0DsCoordinates.col(A);
		vB = GeoClassI.Cell0DsCoordinates.col(B);
		vC = GeoClassI.Cell0DsCoordinates.col(C);
		
		midpoint = ComputePointOnTriangle(3, 1, 1, vA, vB, vC);
		addVertex(GeoClassII, nextVertexId++, midpoint);
		// edgeToFaces
		for (unsigned int edge : GeoClassI.Cell2DsEdges[face])
		{
			edgeToFaces[edge].push_back(face);
		}
		//cout << endl;
		
	}
	
	span<unsigned int> FacesMidpoints(GeoClassII.Cell0DsId.data() + GeoClassI.NumCell0Ds + NumPrincipalEdges,
										GeoClassI.NumCell2Ds);

	
	// connetto i 4 triangoli piccoli per ogni lato principale
	
	unsigned int nextEdgeId = 0;
	unsigned int nextFaceId = 0;
	unsigned int face1;
	unsigned int face2;
	unsigned int faceMidpoint1;
	unsigned int faceMidpoint2;
	unsigned int extrema0;
	unsigned int extrema1;
	unsigned int edgeMidpoint;
	for (int e = 0; e < NumPrincipalEdges; ++e)
	{
		
		extrema0 = GeoClassI.Cell1DsExtrema(0, e);
		extrema1 = GeoClassI.Cell1DsExtrema(1, e);
		face1 = edgeToFaces[e][0];
		face2 = edgeToFaces[e][1];
		faceMidpoint1 = FacesMidpoints[face1];
		faceMidpoint2 = FacesMidpoints[face2];
		edgeMidpoint = EdgesMidpoints[e];
		addFace(GeoClassII, nextEdgeId, nextFaceId, extrema0, edgeMidpoint, faceMidpoint1);
		addFace(GeoClassII, nextEdgeId, nextFaceId, extrema0, edgeMidpoint, faceMidpoint2);
		addFace(GeoClassII, nextEdgeId, nextFaceId, extrema1, edgeMidpoint, faceMidpoint1);
		addFace(GeoClassII, nextEdgeId, nextFaceId, extrema1, edgeMidpoint, faceMidpoint2);
	}
	// connetto i 2 triangoli equilateri a cavallo di ogni lato interno
	for (int e = NumPrincipalEdges; e < GeoClassI.NumCell1Ds; ++e)
	{
		extrema0 = GeoClassI.Cell1DsExtrema(0, e);
		extrema1 = GeoClassI.Cell1DsExtrema(1, e);
		face1 = edgeToFaces[e][0];
		face2 = edgeToFaces[e][1];
		faceMidpoint1 = FacesMidpoints[face1];
		faceMidpoint2 = FacesMidpoints[face2];
		addFace(GeoClassII, nextEdgeId, nextFaceId, extrema0, faceMidpoint1, faceMidpoint2);
		addFace(GeoClassII, nextEdgeId, nextFaceId, extrema1, faceMidpoint1, faceMidpoint2);
		
	}
	NormalizeMatrixColumns(GeoClassII.Cell0DsCoordinates);
	
	return GeoClassII;
}

GeodesicPolyhedron dualize(const GeodesicPolyhedron& poly)
{
    GeodesicPolyhedron dual;
    
    // Il duale ha: vertici = facce originali, facce = vertici originali
    dual.NumCell0Ds = poly.NumCell2Ds;
    dual.NumCell2Ds = poly.NumCell0Ds;
    dual.NumCell1Ds = poly.NumCell1Ds;
    // Inizializza strutture
    InitializeGeodesicStorage(dual);
    
    dual.Cell0DsId = poly.Cell2DsId;
    dual.Cell2DsId = poly.Cell0DsId;
    dual.Cell1DsId = poly.Cell1DsId;
    

    
    // 1. Nuovi vertici = baricentri delle facce originali
    
    
    vector<vector<unsigned int>> edgeToFaces;
	edgeToFaces.resize(poly.NumCell1Ds); // dai lati del geodesico alle facce (corrisponde a lati Goldberg)
	
	int A;
	int B;
	int C;
	
	VectorXd vA;
	VectorXd vB;
	VectorXd vC;
	VectorXd midpoint;
	
	for(unsigned int face : poly.Cell2DsId)
	{
		// punto medio
		A = poly.Cell2DsVertices[face][0];
		B = poly.Cell2DsVertices[face][1];
		C = poly.Cell2DsVertices[face][2];
		
		vA = poly.Cell0DsCoordinates.col(A);
		vB = poly.Cell0DsCoordinates.col(B);
		vC = poly.Cell0DsCoordinates.col(C);
		
		midpoint = ComputePointOnTriangle(3, 1, 1, vA, vB, vC);
		addVertex(dual, face, midpoint);
		// edgeToFaces
		for (unsigned int edge : poly.Cell2DsEdges[face])
		{
			edgeToFaces[edge].push_back(face);
		}
		
	}
    
    // 2. Nuovi lati = connessioni tra facce originali adiacenti
    
    for(unsigned int edge : dual.Cell1DsId) {
	    addEdge(dual, edge, edgeToFaces[edge][0], edgeToFaces[edge][1]);
        }
    
    // 3. Nuove facce = una per ogni vertice originale
    
    for(unsigned int f : dual.Cell2DsId) {
        const vector<unsigned int>& faceEdges = poly.VertextoEdges[f];
        
        unsigned int currentEdge = faceEdges[0];
		unsigned int e1 = dual.Cell1DsExtrema(0, currentEdge);
		unsigned int e2 = dual.Cell1DsExtrema(1, currentEdge);
		
		vector<unsigned int> orderedEdges = { currentEdge };
		vector<unsigned int> orderedVertices = { e1, e2 };
		
		while (true) {
			unsigned int nextStart;
			int nextEdge = findNextEdge(currentEdge, e1, e2, faceEdges, dual.Cell1DsExtrema, nextStart);
		
			if (nextEdge == -1 || nextStart == orderedVertices.front()) break;
		
			orderedEdges.push_back(nextEdge);
			orderedVertices.push_back(
				(dual.Cell1DsExtrema(0, nextEdge) == e2) ? dual.Cell1DsExtrema(1, nextEdge) : dual.Cell1DsExtrema(0, nextEdge)
			);
		
			currentEdge = nextEdge;
			e1 = e2;
			e2 = orderedVertices.back();
		}
		orderedVertices.pop_back();

		dual.Cell2DsVertices[f] = orderedVertices;
		dual.Cell2DsEdges[f] = orderedEdges;
	}
	NormalizeMatrixColumns(dual.Cell0DsCoordinates);
	return dual;
}



int findNextEdge(
    unsigned int currentEdge,
    unsigned int e1,
    unsigned int e2,
    const vector<unsigned int>& faceEdges,
    const MatrixXi& Cell1DsExtrema,
    unsigned int& nextStart)
{
    for (unsigned int edge : faceEdges) {
        if (edge == currentEdge) continue;

        unsigned int a = Cell1DsExtrema(0, edge);
        unsigned int b = Cell1DsExtrema(1, edge);

        if (a == e2 || b == e2) {
            nextStart = (a == e2) ? a : b;
            return edge;
        }
    }

    // Nessun edge trovato (in teoria non dovrebbe accadere se la faccia è chiusa correttamente)
    return -1;
}

 void ComputeShortestPath(GeodesicPolyhedron& mesh, unsigned int source, unsigned int target) {
    const int n = mesh.NumCell0Ds;
    
    // Verifica validità degli indici
    if (source >= n || target >= n) {
        cerr << "Errore: indici di vertice non validi. Source: " << source 
                  << ", Target: " << target << ", NumVertices: " << n << "\n";
        return;
    }
    
    // Inizializzazione delle strutture dati
    vector<double> dist(n, std::numeric_limits<double>::infinity());
    vector<int> pred(n, -1);
    vector<bool> visited(n, false);
    
    // Costruzione della lista di adiacenza per efficienza
    vector<vector<pair<int, double>>> adj(n);
    
    for (unsigned int e = 0; e < mesh.NumCell1Ds; ++e) {
        int v0 = mesh.Cell1DsExtrema(0, e);
        int v1 = mesh.Cell1DsExtrema(1, e);
        
        // Calcola la distanza euclidea tra i vertici
        double edge_length = (mesh.Cell0DsCoordinates.col(v0) - mesh.Cell0DsCoordinates.col(v1)).norm();
        
        // Aggiungi arco in entrambe le direzioni (grafo non diretto)
        adj[v0].emplace_back(v1, edge_length);
        adj[v1].emplace_back(v0, edge_length);
    }
    
    // Algoritmo di Dijkstra
    using P = std::pair<double, int>; // (distanza, nodo)
    priority_queue<P, std::vector<P>, std::greater<P>> pq;
    
    dist[source] = 0.0;
    pq.push({0.0, source});
    
    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        
        if (visited[u]) continue;
        visited[u] = true;
        
        // Se abbiamo raggiunto il target, possiamo terminare
        if (u == target) break;
        
        // Esplora tutti i vicini
        for (const auto& [v, weight] : adj[u]) {
            if (!visited[v] && dist[v] > dist[u] + weight) {
                dist[v] = dist[u] + weight;
                pred[v] = u;
                pq.push({dist[v], v});
            }
        }
    }
    
    // Verifica se esiste un percorso
    if (pred[target] == -1 && target != source) {
        cerr << "Nessun percorso trovato tra " << source << " e " << target << "\n";
        return;
    }
    
    // Inizializza i flag del percorso
    mesh.Cell0DsShortPath.assign(n, false);
    mesh.Cell1DsShortPath.assign(mesh.NumCell1Ds, false);
    
    // Ricostruisce il percorso
    vector<unsigned int> path;
    unsigned int v = target;
    
    while (v != source) {
        path.push_back(v);
        if (pred[v] == -1) break; // Sicurezza aggiuntiva
        v = pred[v];
    }
    path.push_back(source);
    
    // Marca i vertici del percorso
    for (unsigned int vertex : path) {
        mesh.Cell0DsShortPath[vertex] = true;
    }
    
    // Marca gli archi del percorso
    unsigned int edge_count = 0;
    for (int i = path.size() - 1; i > 0; --i) {
        unsigned int u = path[i];
        unsigned int v = path[i-1];
        
        // Trova l'arco tra u e v
        for (unsigned int e = 0; e < mesh.NumCell1Ds; ++e) {
            int a = mesh.Cell1DsExtrema(0, e);
            int b = mesh.Cell1DsExtrema(1, e);
            if ((a == u && b == v) || (a == v && b == u)) {
                mesh.Cell1DsShortPath[e] = true;
                ++edge_count;
                break;
            }
        }
    }
    
    // Output dei risultati
    cout << "Cammino minimo trovato!\n";
    cout << "Numero di lati: " << edge_count << "\n";
    cout << "Lunghezza totale: " << dist[target] << "\n";
    
    // Stampa il percorso
    cout << "Percorso: ";
    for (int i = path.size() - 1; i >= 0; --i) {
        std::cout << path[i];
        if (i > 0) std::cout << " -> ";
    }
    cout << "\n";

}



}



