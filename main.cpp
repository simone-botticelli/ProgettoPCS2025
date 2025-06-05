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


int main() {
	PlatonicType type = PlatonicType::ICOSAHEDRON;
	PlatonicSolid solid(type);
	
	GeodesicPolyhedron geo_1 = PolyhedralLibrary::Build_ClassI_Geodesic(solid, 40);
	GeodesicPolyhedron geo = dualize(geo_1);
	
	ComputeShortestPath(geo, 36, 1000);
	
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





return 0;


}