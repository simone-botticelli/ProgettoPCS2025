#pragma once
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace PolyhedralLibrary;


void ExportCell0Ds(const GeodesicPolyhedron& polyhedron, const string& filePath);
void ExportCell1Ds(const GeodesicPolyhedron& polyhedron, const string& filePath);
void ExportCell2Ds(const GeodesicPolyhedron& polyhedron, const string& filePath);
void ExportCell3Ds(const GeodesicPolyhedron& polyhedron, const string& filePath);
void ExportAllCells(const GeodesicPolyhedron& polyhedron, const string& basePath);