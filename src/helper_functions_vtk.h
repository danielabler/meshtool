#ifndef HELPERFCTS_VTK_H_
#define HELPERFCTS_VTK_H_


#include <iostream>
#include <list>


// VTK
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCellData.h>
#include <vtkIntArray.h>
#include <vtkPolyData.h>
#include <vtkCell.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkPointDataToCellData.h>
#include <vtkDoubleArray.h>
#include <vtkTetra.h>
#include <vtkNIFTIImageReader.h>
#include <vtkNIFTIImageWriter.h>
#include <vtkImageReader2Factory.h>
#include <vtkMetaImageWriter.h>
#include <vtkCellDataToPointData.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkThreshold.h>
#include <vtkPolyData.h>
#include <vtkImageWeightedSum.h>
#include <vtkImageCast.h>

//! @brief hfvtk : Namespace for HELPERFTCS VTK
namespace hfvtk {

    vtkSmartPointer<vtkImageData> loadImage(std::string path_to_image_file);
    vtkSmartPointer<vtkImageData> loadVtkImageData(std::string path_to_file);
    vtkSmartPointer<vtkUnstructuredGrid> loadVtkUnstructuredGrid(std::string path_to_mesh_file);
    vtkSmartPointer<vtkPolyData> extractSurfaceFromVtkDataSet(vtkSmartPointer<vtkDataSet> vtkdataset);
    std::list<int> getBlockIdList(vtkSmartPointer<vtkDataSet> dataset);
    std::list<int> getUniqeIntList(vtkSmartPointer<vtkDataSet> dataset, std::string array_name);
    vtkSmartPointer<vtkImageData> createVtkImageFromVtkPolyDataList(std::list<vtkSmartPointer<vtkPolyData> > vtkPolyDataList, double bounds[6], double spacing[3]);
    vtkSmartPointer<vtkImageData> createVtkImageFromVtkPolyData(vtkSmartPointer<vtkPolyData> vtkPolyData, double bounds[6], double spacing[3]);
    vtkSmartPointer<vtkImageData> createVtkImageFromVtkUnstructuredGrid(vtkSmartPointer<vtkUnstructuredGrid> vtk_ugrid, double spacing[3]);
    void writeVtkImageData(vtkSmartPointer<vtkImageData> image_data, std::string path_to_output_file);
    void writeVtkUnstructuredGrid(vtkSmartPointer<vtkUnstructuredGrid> , std::string );
    void writeVtkData(vtkSmartPointer<vtkDataSet>, std::string );

}
#endif // HELPERFCTS_VTK_H_