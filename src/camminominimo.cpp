#include <vector>
#include <queue>
#include <limits>
#include <iostream>
#include <cmath>
#include <unordered_map>
#include "PolyhedralMesh.hpp"

using namespace std;

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
        int v0 = mesh.Cell1DsExtrema(e, 0);
        int v1 = mesh.Cell1DsExtrema(e, 1);
        }
        
        // Calcola la distanza euclidea tra i vertici
        double dx = mesh.Cell0DsCoordinates(v0, 0) - mesh.Cell0DsCoordinates(v1, 0);
        double dy = mesh.Cell0DsCoordinates(v0, 1) - mesh.Cell0DsCoordinates(v1, 1);
        double dz = mesh.Cell0DsCoordinates(v0, 2) - mesh.Cell0DsCoordinates(v1, 2);
        double edge_length = sqrt(dx*dx + dy*dy + dz*dz);
        
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
            int a = mesh.Cell1DsExtrema(e, 0);
            int b = mesh.Cell1DsExtrema(e, 1);
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