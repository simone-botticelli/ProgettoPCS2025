#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>

#include "PolyhedralMesh.hpp"
#include "UCDUtilities.hpp" 
#include "ExportTxt.hpp"


using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;
using namespace Gedim;

int main(int argc, char** argv) {
    // Verifica che ci siano esattamente 5 argomenti (programma + p, q, b1, b2)
    if (argc != 5 && argc != 7) {
        cerr << "Uso corretto: ./ProgettoPCS <p> <q> <b> <c> oppure ./ProgettoPCS <p> <q> <b> <c> <id_start> <id_end> \n";
        return 1;
    }

int p, q, b, c;

	try {
    p = stoi(argv[1]);
    q = stoi(argv[2]);
    b = stoi(argv[3]);
    c = stoi(argv[4]);
	} catch (const exception& e) {
    	cerr << "Errore: uno degli argomenti <p> <q> <b> <c> non è un intero valido.\n";
    	return 1;
	}
    
    bool to_dualize = false;
    
    if ((p != 3 && q != 3) || p < 3 || p > 5 || q < 3 || q > 5) {
	    cerr << "Errore: almeno uno tra p o q deve essere 3, l'altro deve essere compreso tra 3 e 5.\n";
	    return 1;
	}
	    
	if ((b != 0 && c != 0 && b != c) || (b == 0 && c == 0)) {
		cerr << "Errore: esattamente uno tra b e c deve essere diverso da 0, oppure b e c devono essere uguali.\n";
		return 1;
	}
	 
    if (p != 3 && q == 3) {
	    to_dualize = true;
	    q = p;
	}
    
    PlatonicType type;
    
	switch(q) {
		case 3:
			type = PlatonicType::TETRAHEDRON;
			break;
		case 4:
			type = PlatonicType::OCTAHEDRON;
			break;
		case 5:
			type = PlatonicType::ICOSAHEDRON;
			break;
    }

	PlatonicSolid solid(type);
	GeodesicPolyhedron geo; 
	
	if (b != c) {
	geo = PolyhedralLibrary::Build_ClassI_Geodesic(solid, b+c);
	cout << "E' stato selezionato il poliedro di classe I.\n";
	}
	else {
	geo = PolyhedralLibrary::Build_ClassII_Geodesic(solid, b);
	cout << "E' stato selezionato il poliedro di classe II.\n";
	}
	
	if (to_dualize) {
		geo = dualize(geo);
		cout << "Ho costruito il poliedro duale come richiesto.\n";
	}
	
if (argc == 7) {
    try {
        int start = stoi(argv[5]);
        int end   = stoi(argv[6]);
        ComputeShortestPath(geo, start, end);
    } catch (const exception& e) {
        cerr << "Argomenti <start> o <end> non validi: devono essere interi.\n";
    }
}
	
	UCDUtilities utilities;
    {
	    vector<UCDProperty<double>> cell0Ds_properties(1);
	
		
        cell0Ds_properties[0].Label = "ShortPath";
        cell0Ds_properties[0].UnitLabel = "-";
        cell0Ds_properties[0].NumComponents = 1;

        vector<double> cell0Ds_marker(geo.NumCell0Ds, 0.0);

        for(int i = 0; i < geo.NumCell0Ds; i++)
        {
                bool m = geo.Cell0DsShortPath[i];
                cell0Ds_marker[i] = double(m);
                
        }

        cell0Ds_properties[0].Data = cell0Ds_marker.data();
                
        utilities.ExportPoints("./Cell0Ds.inp",
                               geo.Cell0DsCoordinates,
                               cell0Ds_properties);
    }
    
    {
	    
	    vector<UCDProperty<double>> cell1Ds_properties(1);
	
		
        cell1Ds_properties[0].Label = "ShortPath";
        cell1Ds_properties[0].UnitLabel = "-";
        cell1Ds_properties[0].NumComponents = 1;

        vector<double> cell1Ds_marker(geo.NumCell1Ds, 0.0);

        for(int i = 0; i < geo.NumCell1Ds; i++)
        {
                bool m = geo.Cell1DsShortPath[i];
                cell1Ds_marker[i] = double(m);
                
        }

        cell1Ds_properties[0].Data = cell1Ds_marker.data();
        utilities.ExportSegments("./Cell1Ds.inp",
                                 geo.Cell0DsCoordinates,
                                 geo.Cell1DsExtrema,
                                 {},
                                 cell1Ds_properties);
    }
ExportAllCells(geo, "./");

cout << "Poliedro geodetico generato con successo!\n";

return 0;


}