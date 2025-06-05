#include <iostream>
#include <cassert>
#include <queue>
#include <limits>
#include <iostream>
#include <cmath>

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
};

GeodesicCounts SetGeodesicCounts_ClassI(const PlatonicSolid& solid, const unsigned int& n)
{
	GeodesicCounts counts;
	unsigned int T = n*n;
	
	switch(solid.Type) 
	{
		case PlatonicType::TETRAHEDRON:
			counts.V = 2*T + 2;
			counts.E = 6*T;
			counts.F = 4*T;
				return counts;
			
		case PlatonicType::OCTAHEDRON:
			counts.V = 4*T + 2;
			counts.E = 12*T;
			counts.F = 8*T;
				return counts;
		
		case PlatonicType::ICOSAHEDRON:
			counts.V = 10*T + 2;
			counts.E = 30*T;
			counts.F = 20*T;
				return counts;
			
		default:
			cerr << "PlatonicType not supported" << endl;
			throw std::runtime_error("Unsupported PlatonicType in SetGeodesicCounts_ClassI");
	}
}

void addVertex(GeodesicPolyhedron& geodesic, unsigned int vertexId, const VectorXd& vertexCoordinates)
{
	assert(vertexId < geodesic.Cell0DsCoordinates.cols() && "vertexId fuori range");
	geodesic.Cell0DsId.push_back(vertexId);
	geodesic.Cell0DsCoordinates.col(vertexId) = vertexCoordinates;
}

void addEdge(GeodesicPolyhedron& geodesic, unsigned int edgeId, unsigned int originId, unsigned int endId)
{
	assert(edgeId < geodesic.Cell1DsExtrema.cols() && "edgeId fuori range");
	geodesic.Cell1DsId.push_back(edgeId);
	geodesic.Cell1DsExtrema(0, edgeId) = originId;
	geodesic.Cell1DsExtrema(1, edgeId) = endId;
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
            geodesic.Cell1DsId.push_back(newId);
            geodesic.Cell1DsExtrema(0, newId) = originId;
			geodesic.Cell1DsExtrema(1, newId) = endId;
			return newId;
            }
       
}

// Overload per triangoli
void addFace(
    GeodesicPolyhedron& geodesic,
    unsigned int& nextEdgeId,
    unsigned int faceId,
    unsigned int v0, unsigned int v1, unsigned int v2)
{
    addFace(geodesic,nextEdgeId, faceId, vector<unsigned int>{v0,v1,v2});
}

void addFace( 
    GeodesicPolyhedron& geodesic,
    unsigned int& nextEdgeId,
    unsigned int faceId,
    const vector<unsigned int>& facesVertices
)
{
    // controlli...
    int N = facesVertices.size();
    vector<unsigned int> facesEdges(N);
    for(int i = 0; i < N; i++){
	    unsigned int edgeId = GetorAddEdge(geodesic, nextEdgeId, facesVertices[i], facesVertices[(i + 1)%N]);
	    facesEdges.push_back(edgeId);
	    
    }
    geodesic.Cell2DsId.push_back(faceId);
    geodesic.Cell2DsVertices.push_back(facesVertices);
    geodesic.Cell2DsEdges.push_back(facesEdges);
}


GeodesicPolyhedron Build_ClassI_Geodesic(const PlatonicSolid& solid, const unsigned int n)
{
	GeodesicPolyhedron geodesic;

	// 1. Inizializzo NumCell0Ds, NumCell1Ds, NumCell2Ds
	
		GeodesicCounts counts = SetGeodesicCounts_ClassI(solid, n);
		geodesic.NumCell0Ds = counts.V;
		geodesic.NumCell1Ds = counts.E;
		geodesic.NumCell2Ds = counts.F;
		// cout << "NumCell0Ds: " << geodesic.NumCell0Ds << endl;
		// cout << "NumCell1Ds: " << geodesic.NumCell1Ds << endl;
		// cout << "NumCell2Ds: " << geodesic.NumCell2Ds << endl;
	
	// 2. Inizializzo le dimensioni dei vettori Cell0DsId, Cell1DsId, Cell2DsId, Cell0DsShortPath, Cell1DsShortPath
		
		geodesic.Cell0DsId.reserve(geodesic.NumCell0Ds);
		geodesic.Cell1DsId.reserve(geodesic.NumCell1Ds);
		geodesic.Cell2DsId.reserve(geodesic.NumCell2Ds);
		geodesic.Cell0DsShortPath.reserve(geodesic.NumCell0Ds);
		geodesic.Cell1DsShortPath.reserve(geodesic.NumCell1Ds);
		
	// 3. Inizializzo le dimensioni delle matrici Cell0DsCoordinates, Cell1DsExtrema, Cell2DsVertices, Cell2DsEdges

		// Matrice coordinate punti 0D: 3 x NumCell0Ds (x,y,z)
		geodesic.Cell0DsCoordinates.resize(3, geodesic.NumCell0Ds);
		
		// Matrice estremi 1D: 2 x NumCell1Ds (due vertici per ogni cella 1D)
		geodesic.Cell1DsExtrema.resize(2, geodesic.NumCell1Ds);
		
		// Ridimensiono vettore esterno con NumCell2Ds elementi
		geodesic.Cell2DsVertices.reserve(geodesic.NumCell2Ds);
		geodesic.Cell2DsEdges.reserve(geodesic.NumCell2Ds);

	// 4. Salvo i dati dei vertici del solido platonico (id e coordinate) in  Cell0DsCoordinates e in Cell0DsId
		
		geodesic.Cell0DsId = solid.VerticesId;  // copia vettore di ID
		geodesic.Cell0DsCoordinates.block(0, 0, solid.VerticesCoordinates.rows(), solid.VerticesCoordinates.cols()) = solid.VerticesCoordinates;
	
	// se n = 1 mi fermo qui.
	
	if (n == 1)
	{
		// copiare struttura solido platonico
		geodesic.Cell1DsId = solid.EdgesId;
		geodesic.Cell2DsId = solid.FacesId;
		geodesic.Cell1DsExtrema = solid.EdgesExtrema;
		geodesic.Cell2DsVertices = solid.FacesVertices;
		geodesic.Cell2DsEdges = solid.FacesEdges;
		geodesic.Cell3DsId = 0;
		
		return geodesic;
	}
	
	// 5.  SUDDIVISIONE LATI PRINCIPALI:  
				  
		// per ogni lato in solid.edges lo suddivido (in modo di partire dal vertice con id più basso)
		
		unsigned int nextVertexId = solid.NumVertices; // ho già inizializzato i vertici principali
		unsigned int oldVertexId = nextVertexId;
		unsigned int nextEdgeId = 0;
		unsigned int oldEdgeId = nextEdgeId;
		
		for(int e = 0; e < solid.NumEdges; e++){

			int v1 = solid.EdgesExtrema(0, e);
			int v2   = solid.EdgesExtrema(1, e);
			
			int v_start = min(v1,v2);
			int v_end = max(v1,v2);

			// ottengo le coordinate dei vertici estremi
			VectorXd start = solid.VerticesCoordinates.col(v_start);
			VectorXd end   = solid.VerticesCoordinates.col(v_end);
			
			// 5.1 Gestisco separatamente il primo lato piccolo
				
				// creo le coordinate del primo vertice
				double t = double(1) / n;
				VectorXd new_point = (1-t)*start + t*end;
				
				// salvo il primo vertice
				addVertex(geodesic, nextVertexId, new_point);
				
				// salvo il primo lato piccolo 
				GetorAddEdge(geodesic, nextEdgeId, v_start, nextVertexId);
				
				oldVertexId = nextVertexId++; // assegna e poi aumenta

			// 5.2 Creo gli altri vertici e lati sul lato principale
			
				for (int k = 2; k < n; ++k)  // k va da 2 a n-1
				{
					
					// creo le coordinate
					double t = double(k) / double(n);
					VectorXd new_point = (1-t)*start + t*end;
					
					// salvo il vertice
					addVertex(geodesic, nextVertexId, new_point);
					
					// salvo gli estremi dei lati in Cell1DsExtrema
					GetorAddEdge(geodesic, nextEdgeId, oldVertexId, nextVertexId);
					
					oldVertexId = nextVertexId++;
					
				}
			
			// 5.3 Collego l'ultimo lato piccolo
				GetorAddEdge(geodesic, nextEdgeId, oldVertexId, v_end);
			}
		
			// Usando Eigen Map non creo copie ma un alias che mi permette di lavorare su una parte del vector originale
		
			// 5.4 Con Eigen Map mi creo una matrice che nell'i-esima colonna ha gli id dei vertici interni al lato i-esimo del solido
				
				
				// Controllo di sicurezza (poi lo leviamo)
				assert(geodesic.Cell0DsId.size() >= solid.NumVertices + solid.NumEdges * (n - 1));
				
				// Mappo direttamente in una matrice Eigen senza copiare
				const Map<Matrix<unsigned int, Dynamic, Dynamic, ColMajor>> 
					internalVerticesMatrix(geodesic.Cell0DsId.data() + solid.NumVertices, // punto di partenza nel vector
										   n - 1, // elementi per ogni colonna 
										   solid.NumEdges); // num colonne
		
			// 5.5 Con Eigen Map mi creo una matrice che nell'i-esima colonna ha gli id dei lati interni al lato i-esimo del solido
			
				const Map<Matrix<unsigned int, Dynamic, Dynamic, ColMajor>> 
					internalEdgesMatrix(geodesic.Cell1DsId.data(), 
										n, 
										solid.NumEdges);
			
			// 5.6 Controlli ...

	// 6. TRIANGOLAZIONE
		unsigned int nextFaceId = 0;
		int i = 1;
		int j = 1;
		double alpha = double(n - i - j) / n;
		double beta  = double(i) / n;
		double gamma = double(j) / n;
		
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
				addFace(geodesic, nextEdgeId, nextFaceId++, 
							A, AB[0], AC[0]);
				addFace(geodesic, nextEdgeId, nextFaceId++, 
							B, AB[0], BC[0]);
				addFace(geodesic, nextEdgeId, nextFaceId++, 
							C, AC[0], BC[0]);
				addFace(geodesic, nextEdgeId, nextFaceId++, 
							BC[0], AB[0], AC[0]);
							
				continue;
				}
							
			// caso n > 2 
				
				// passo 1: 
			
					// salvo la prima faccia (ATTENZIONE all ordine di id_vertici e  id_lati)
					addFace(geodesic, nextEdgeId, nextFaceId++, 
							A, AB[0], AC[0]);
				
				// passo 2: creo il primo punto interno
					
					i = 1;
					j = 1;
					alpha = double(n - i - j) / n;
					beta  = double(i) / n;
					gamma = double(j) / n;
					
					VectorXd new_point = alpha * vA + beta * vB + gamma * vC;
					addVertex(geodesic, nextVertexId, new_point);
					
					
					// 2.1 salvo la prima faccia
						addFace(geodesic, nextEdgeId, nextFaceId++,
								AB[0], AB[1], nextVertexId);

					// 2.2 salvo la seconda faccia
						addFace(geodesic, nextEdgeId, nextFaceId++,
								nextVertexId, AC[0], AB[0]);
					
					// 2.3 salvo la terza faccia 
						addFace(geodesic, nextEdgeId, nextFaceId++,
								nextVertexId, AC[1], AC[0]);
							
					nextVertexId++;
										 
				// passo 3:  
				
				// per ogni riga (tranne l'ultima)
				for (int s = 3; s < n; s++)
				{
				// 3.1 creo il primo punto interno (j = 1)
					// Coordinate baricentriche
					j = 1;
					i = s - j;
					alpha = double(n - i - j) / n;
					beta  = double(i) / n;
					gamma = double(j) / n;
					

					new_point = alpha * vA + beta * vB + gamma * vC;
					addVertex(geodesic, nextVertexId, new_point);
						
						// salvo faccia a sx 
						addFace(geodesic, nextEdgeId, nextFaceId++,
								AB[i - 1], AB[i], nextVertexId);
						
						// salvo faccia sotto
						addFace(geodesic, nextEdgeId, nextFaceId++,
								AB[i - 1], nextVertexId, nextVertexId - i + 1);

					nextVertexId++;
						
					// 3.2 creo il secondo vertice interno
					
					j = 2;
					i = s - j;
					alpha = double(n - i - j) / n;
					beta  = double(i) / n;
					gamma = double(j) / n;
					new_point = alpha * vA + beta * vB + gamma * vC;
					addVertex(geodesic, nextVertexId, new_point); 
						
						
						// creo la faccia a sx 
						addFace(geodesic, nextEdgeId, nextFaceId++,
								nextVertexId - 1, nextVertexId, nextVertexId - s + 1);
						
					nextVertexId++;
					
					// 3.3 per ogni colonna (da 2 alla penultima)
					for (j = 3; j < s; j++)
					{
						
						i = s - j;

						// creo la faccia sotto 
						addFace(geodesic, nextEdgeId, nextFaceId++,
								nextVertexId - 1, nextVertexId - s + 1, nextVertexId - s);
						

						// aggiungo il vertice della colonna j + 1 
						alpha = double(n - i - j) / n;
						beta  = double(i) / n;
						gamma = double(j) / n;
						VectorXd new_point = alpha * vA + beta * vB + gamma * vC;
						addVertex(geodesic, nextVertexId, new_point);

						// creo la faccia alla sua sx
						addFace(geodesic, nextEdgeId, nextFaceId++,
								nextVertexId - s + 1, nextVertexId - 1, nextVertexId);
						
						nextVertexId++;
					}
						
					// 3.4 ultima colonna
						
						// collego l'ultimo vertice ad AC[j - 1] e creo la faccia sotto
						
						addFace(geodesic, nextEdgeId, nextFaceId++,
								AC[s-2], nextVertexId - s , nextVertexId - 1);
						
						// collego l'ultimo vertice ad AC[j] e creo la faccia a dx
						
						addFace(geodesic, nextEdgeId, nextFaceId++,
								AC[s-1], AC[s-2], nextVertexId - 1);
						
				}
				// passo 4: ultima riga
					 
					// 4.1 collego le prime tre facce 
					// collego AB[n-2]-BC[0] e creo la prima faccia
						
						
						addFace(geodesic, nextEdgeId, nextFaceId++,
								AB[n-2], B, BC[0]);
						
			
						// collego BC[0] al vertice sotto e creo la seconda faccia
						
						addFace(geodesic, nextEdgeId, nextFaceId++,
								AB[n-2], BC[0], nextVertexId - n + 2);
						
						// collego BC[1] al vertice sotto a sx e creo la terza faccia
						
						addFace(geodesic, nextEdgeId, nextFaceId++,
								nextVertexId - n + 2, BC[0], BC[1]);
								
					// 4.2 per ogni vertice sul lato BC (dal secondo al penultimo[1:n-3]) 
					for(int k = 1; k < n - 2; k++)
					{
						// lo collego al punto sotto a dx e creo una faccia
						int vertex_under = nextVertexId - n + 2 + k;
						addFace(geodesic, nextEdgeId, nextFaceId++,
								vertex_under - 1, BC[k], vertex_under);

						// collego il vertice successivo a quello sotto e creo un altra faccia
						
						addFace(geodesic, nextEdgeId, nextFaceId++,
								BC[k], BC[k + 1], vertex_under);
						
					}
						
					// 4.3 collego il lato BC[n-2]-AC[n-2] e salvo le ultime 2 facce
						
						addFace(geodesic, nextEdgeId, nextFaceId++,
								BC[n-2], AC[n-2], nextVertexId - 1);
						addFace(geodesic, nextEdgeId, nextFaceId++,
								BC[n-2], C, AC[n-2]);
	
		}  // fine ciclo sulle facce 
		  
	// Normalizzazione
	MatrixXd& M = geodesic.Cell0DsCoordinates;  // MatrixXd M(n, m) già inizializzata
		// Normalizza ogni colonna in-place
		for (int j = 0; j < M.cols(); ++j) {
			double norm = M.col(j).norm();
			if (norm > 1e-12) {
				M.col(j) /= norm;
			}
		}

	return geodesic;
	
} // fine build_classI_Geodesic
	

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

GeodesicPolyhedron dualize(const GeodesicPolyhedron& poly)
{
    GeodesicPolyhedron dual;
    
    // Il duale ha: vertici = facce originali, facce = vertici originali
    dual.NumCell0Ds = poly.NumCell2Ds;
    dual.NumCell2Ds = poly.NumCell0Ds;
    dual.NumCell1Ds = poly.NumCell1Ds;
    
    // Inizializza strutture
    dual.Cell0DsId.resize(dual.NumCell0Ds);
    dual.Cell2DsId.resize(dual.NumCell2Ds);
    dual.Cell1DsId.resize(dual.NumCell1Ds);
    for(unsigned int i = 0; i < dual.NumCell0Ds; i++) dual.Cell0DsId[i] = i;
    for(unsigned int i = 0; i < dual.NumCell2Ds; i++) dual.Cell2DsId[i] = i;
    for(unsigned int i = 0; i < dual.NumCell1Ds; i++) dual.Cell1DsId[i] = i;
    
    dual.Cell0DsShortPath.resize(dual.NumCell0Ds, false);
    dual.Cell1DsShortPath.resize(dual.NumCell1Ds, false);
    dual.Cell3DsId = 0;
    
    // 1. Nuovi vertici = baricentri delle facce originali
    dual.Cell0DsCoordinates.resize(3, dual.NumCell0Ds);
    
    for(unsigned int f = 0; f < poly.NumCell2Ds; f++) {
        VectorXd centroid = VectorXd::Zero(3);
        
        // Calcola baricentro della faccia f
        for(unsigned int v : poly.Cell2DsVertices[f]) {
            centroid += poly.Cell0DsCoordinates.col(v);
        }
        centroid /= poly.Cell2DsVertices[f].size();
        
        // Proietta sulla sfera unitaria
        dual.Cell0DsCoordinates.col(f) = centroid.normalized();
    }
    
    // 2. Nuovi lati = connessioni tra facce originali adiacenti
    dual.Cell1DsExtrema.resize(2, dual.NumCell1Ds);
    
    for(unsigned int e = 0; e < poly.NumCell1Ds; e++) {
        // Trova le 2 facce che condividono il lato e
        vector<unsigned int> faces_with_edge;
        
        for(unsigned int f = 0; f < poly.NumCell2Ds; f++) {
            for(unsigned int edge : poly.Cell2DsEdges[f]) {
                if(edge == e) {
                    faces_with_edge.push_back(f);
                    break;
                }
            }
        }
        
        // Connetti i baricentri di queste 2 facce
        if(faces_with_edge.size() == 2) {
            dual.Cell1DsExtrema(0, e) = faces_with_edge[0];
            dual.Cell1DsExtrema(1, e) = faces_with_edge[1];
        }
    }
    
    // 3. Nuove facce = una per ogni vertice originale
    dual.Cell2DsVertices.resize(dual.NumCell2Ds);
    dual.Cell2DsEdges.resize(dual.NumCell2Ds);
    
    for(unsigned int v = 0; v < poly.NumCell0Ds; v++) {
        // Trova tutte le facce che contengono il vertice v
        vector<unsigned int> faces_around_vertex;
        
        for(unsigned int f = 0; f < poly.NumCell2Ds; f++) {
            for(unsigned int vertex : poly.Cell2DsVertices[f]) {
                if(vertex == v) {
                    faces_around_vertex.push_back(f);
                    break;
                }
            }
        }
        
        dual.Cell2DsVertices[v] = faces_around_vertex;
        
        // Trova i lati della nuova faccia
        vector<unsigned int> face_edges;
        for(unsigned int e = 0; e < dual.NumCell1Ds; e++) {
            unsigned int v1 = dual.Cell1DsExtrema(0, e);
            unsigned int v2 = dual.Cell1DsExtrema(1, e);
            
            // Se entrambi i vertici del lato sono nella faccia
            bool has_v1 = false, has_v2 = false;
            for(unsigned int fv : faces_around_vertex) {
                if(fv == v1) has_v1 = true;
                if(fv == v2) has_v2 = true;
            }
            if(has_v1 && has_v2) {
                face_edges.push_back(e);
            }
        }
        dual.Cell2DsEdges[v] = face_edges;
    }
	return dual;
}

}



