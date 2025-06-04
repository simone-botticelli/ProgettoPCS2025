#include <iostream>
#include <cstdlib>
#include "PolyhedralMesh.hpp"
#include "UCDUtilities.hpp" 
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;
using namespace Gedim;




// Funzione per esportare le celle 0D (vertici) in formato .txt
void ExportCell0Ds(const GeodesicPolyhedron& polyhedron, const string& filePath)
{
    ofstream file(filePath);
    if (!file.is_open()) {
        throw runtime_error("Cannot open file: " + filePath);
    }

    file << "# Cell0Ds - Vertices" << endl;
    file << "# Id X Y Z" << endl;
    file << polyhedron.NumCell0Ds << endl;

    file << fixed << setprecision(16);
    for (unsigned int i = 0; i < polyhedron.NumCell0Ds; i++) {
        file << polyhedron.Cell0DsId[i] << " "
             << polyhedron.Cell0DsCoordinates(0, i) << " "
             << polyhedron.Cell0DsCoordinates(1, i) << " "
             << polyhedron.Cell0DsCoordinates(2, i) << endl;
    }

    file.close();
}

// Funzione per esportare le celle 1D (lati) in formato .txt
void ExportCell1Ds(const GeodesicPolyhedron& polyhedron, const string& filePath)
{
    ofstream file(filePath);
    if (!file.is_open()) {
        throw runtime_error("Cannot open file: " + filePath);
    }

    file << "# Cell1Ds - Edges" << endl;
    file << "# Id Origin End" << endl;
    file << polyhedron.NumCell1Ds << endl;

    for (unsigned int i = 0; i < polyhedron.NumCell1Ds; i++) {
        file << polyhedron.Cell1DsId[i] << " "
             << polyhedron.Cell1DsExtrema(0, i) << " "
             << polyhedron.Cell1DsExtrema(1, i) << endl;
    }

    file.close();
}

// Funzione per esportare le celle 2D (facce) in formato .txt
void ExportCell2Ds(const GeodesicPolyhedron& polyhedron, const string& filePath)
{
    ofstream file(filePath);
    if (!file.is_open()) {
        throw runtime_error("Cannot open file: " + filePath);
    }

    file << "# Cell2Ds - Faces" << endl;
    file << "# Id NumVertices NumEdges Vertices Edges" << endl;
    file << polyhedron.NumCell2Ds << endl;

    for (unsigned int i = 0; i < polyhedron.NumCell2Ds; i++) {
        file << polyhedron.Cell2DsId[i] << " "
             << polyhedron.Cell2DsVertices[i].size() << " "
             << polyhedron.Cell2DsEdges[i].size();

        // Esporta i vertici della faccia
        for (unsigned int v : polyhedron.Cell2DsVertices[i]) {
            file << " " << v;
        }

        // Esporta i lati della faccia
        for (unsigned int e : polyhedron.Cell2DsEdges[i]) {
            file << " " << e;
        }

        file << endl;
    }

    file.close();
}

// Funzione per esportare le celle 3D (poliedro) in formato .txt
void ExportCell3Ds(const GeodesicPolyhedron& polyhedron, const string& filePath)
{
    ofstream file(filePath);
    if (!file.is_open()) {
        throw runtime_error("Cannot open file: " + filePath);
    }

    file << "# Cell3Ds - Polyhedron" << endl;
    file << "# Id NumVertices NumEdges NumFaces" << endl;
    file << "0" << endl;  // Un solo poliedro

    file << polyhedron.Cell3DsId << " "
         << polyhedron.NumCell0Ds << " "
         << polyhedron.NumCell1Ds << " "
         << polyhedron.NumCell2Ds << endl;

    file.close();
}

// Funzione per esportare tutti i file .txt richiesti
void ExportAllCells(const GeodesicPolyhedron& polyhedron, const string& basePath)
{
    ExportCell0Ds(polyhedron, basePath + "Cell0Ds.txt");
    ExportCell1Ds(polyhedron, basePath + "Cell1Ds.txt");
    ExportCell2Ds(polyhedron, basePath + "Cell2Ds.txt");
    ExportCell3Ds(polyhedron, basePath + "Cell3Ds.txt");
}




int main() {
	PlatonicType type = PlatonicType::TETRAHEDRON;
	PlatonicSolid solid(type);
	
	GeodesicPolyhedron geo = PolyhedralLibrary::Build_ClassI_Geodesic(solid, 7);
	
	UCDUtilities utilities;
    {
        utilities.ExportPoints("./Cell0Ds.inp",
                               geo.Cell0DsCoordinates);
    }
    
    {
        utilities.ExportSegments("./Cell1Ds.inp",
                                 geo.Cell0DsCoordinates,
                                 geo.Cell1DsExtrema);
    }
ExportAllCells(geo, "./");


return 0;


}