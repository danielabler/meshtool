// file convertmplex_3_subdomain_in_triangulation_3_to_vtk.h

/*!
 * 2015 Daniel Abler
 * 
 * This file uses parts of Andrew Slaughter's 
 * Complex_3_subdomain_in_triangulation_3_to_vtk.h
 * http://aeslaughter.github.io/postdoc/doc/index.html 
 */

#ifndef CGALTOVTK_H_
#define CGALTOVTK_H_

// Standard includes
#include <map>

// VTK includes
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellArray.h>
#include <vtkType.h>

#include <vtkCellData.h>
#include <vtkIntArray.h>
#include <vtkSmartPointer.h>


template <typename C3T3>
vtkSmartPointer<vtkUnstructuredGrid>
c3t3_subdomain_to_vtk_unstructured_grid(const C3T3& c3t3)
{
  
	typedef typename C3T3::Triangulation Triangulation;
	typedef typename Triangulation::Vertex_handle Vertex_handle;

	const Triangulation& tr = c3t3.triangulation();
	
	vtkSmartPointer<vtkPoints> const vtk_points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> const vtk_facets = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkCellArray> const vtk_cells = vtkSmartPointer<vtkCellArray>::New();

	// Create an array for storing the subdomain index
	vtkSmartPointer<vtkIntArray> blockIds = vtkSmartPointer<vtkIntArray>::New();
	blockIds->SetName("ElementBlockIds");

	vtk_points->Allocate(c3t3.triangulation().number_of_vertices());
	vtk_facets->Allocate(c3t3.number_of_facets_in_complex());
	vtk_cells->Allocate(c3t3.number_of_cells_in_complex());

	std::map<Vertex_handle, vtkIdType> V;
	vtkIdType inum = 0;

	for(typename Triangulation::Finite_vertices_iterator 
		vit = tr.finite_vertices_begin(),
		end = tr.finite_vertices_end();
	  vit != end;
	  ++vit)
	{
	typedef typename Triangulation::Point Point;
	const Point& p = vit->point();
	vtk_points->InsertNextPoint(CGAL::to_double(p.x()),
				    CGAL::to_double(p.y()),
				    CGAL::to_double(p.z()));
	V[vit] = inum++;
	}
	
	for(typename C3T3::Facets_in_complex_iterator 
		fit = c3t3.facets_in_complex_begin(),
		end = c3t3.facets_in_complex_end();
	  fit != end; ++fit) 
	{
	vtkIdType cell[3];
	int j=0;
	for (int i = 0; i < 4; ++i)
	  if (i != fit->second)
		cell[j++] =  V[(*fit).first->vertex(i)];
	CGAL_assertion(j==3);
	vtk_facets->InsertNextCell(3, cell);
	}
	for(typename C3T3::Cells_in_complex_iterator 
		cit = c3t3.cells_in_complex_begin(),
		end = c3t3.cells_in_complex_end();
	  cit != end; ++cit) 
	{
	  
	// Add the subdomain index to the vtk array
	blockIds->InsertNextValue( (int)c3t3.subdomain_index(cit) );	  
	  
	vtkIdType cell[4];
	for (int i = 0; i < 4; ++i)
	  cell[i] =  V[cit->vertex(i)];
	vtk_cells->InsertNextCell(4, cell);
	}

	//if(!grid) {
	vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	//}
	
	grid->SetPoints(vtk_points);
	grid->SetCells(VTK_TRIANGLE, vtk_facets);
	grid->SetCells(VTK_TETRA, vtk_cells);
	
	// Add the subdomain index to the vtk grid
	grid->GetCellData()->AddArray(blockIds);
	
	return grid;
        

 }

#endif  // CGALTOVTK_H_


