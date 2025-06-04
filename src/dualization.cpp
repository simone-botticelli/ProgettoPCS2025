#include <vector>
#include <Eigen/Dense>
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace Eigen;

GeodesicPolyhedron dualize(const GeodesicPolyhedron& poly) {
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
        Vector3d centroid = Vector3d::Zero();
        
        // Calcola baricentro della faccia f
        for(unsigned int v : poly.Cell2DsVertices[f]) {
            centroid += poly.Cell0DsCoordinates.col(v);
        }
        centroid /= poly.Cell2DsVertices[f].size();
        
        // Proietta sulla sfera unitaria
        dual.Cell0DsCoordinates.col(f) = centroid.normalized();
    }
    
    // 2. Nuovi lati = connessioni tra facce originali adiacenti
    dual.Cell1DsExtrema.resize(dual.NumCell1Ds, 2);
    
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
            dual.Cell1DsExtrema(e, 0) = faces_with_edge[0];
            dual.Cell1DsExtrema(e, 1) = faces_with_edge[1];
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
            unsigned int v1 = dual.Cell1DsExtrema(e, 0);
            unsigned int v2 = dual.Cell1DsExtrema(e, 1);
            
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