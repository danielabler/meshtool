/**
 * @file   divhelpers.cpp
 * @author Daniel Abler <daniel.abler@istb.unibe.ch>
 * @date   July 2015
 * @version 1.0
 * @brief  Helper Functions
 */


#include "helper_functions_vtk.h"
// LOGGING
#include "easylogging++.h"
#include "helper_functions.h"

namespace hfvtk {


/**
* Function to load image and to store to internal variable
*/
vtkSmartPointer<vtkImageData> loadImage(std::string path_to_image_file)
{
    // initialise images
    vtkSmartPointer<vtkImageData> vtk_image = vtkSmartPointer<vtkImageData>::New();
    std::string inr_suffix = "inr";
    std::string inr_suffix_gz = "inr.gz";
    std::string vti_suffix = ".vti";
    std::string nii_suffix = ".nii";

    if (hf::existsFile(path_to_image_file))
    {
        if ( hf::hasSuffix(path_to_image_file, inr_suffix) or hf::hasSuffix(path_to_image_file, inr_suffix_gz))
        {
            LOG(INFO) << "Input is INRIA file -- currently not supported";
        }
        else
        {
            if ( hf::hasSuffix(path_to_image_file, vti_suffix) )
            {
                vtk_image = loadVtkImageData(path_to_image_file);
            }

            else if ( hf::hasSuffix(path_to_image_file, nii_suffix) )
            {
                LOG(DEBUG)<< "Attempting to load image '" << path_to_image_file <<"' as nifti with vtkNIFTIImageReader";
                vtkSmartPointer<vtkNIFTIImageReader> imageReader =  vtkSmartPointer<vtkNIFTIImageReader>::New();
                imageReader->SetFileName(path_to_image_file.c_str());
                imageReader->Update();
                vtk_image = imageReader->GetOutput();
            }
            else
            {
                LOG(DEBUG)<< "Attempting to load image '" << path_to_image_file <<"' as other vtk format with vtkImageReader2Factory";
                vtkSmartPointer<vtkImageReader2Factory> readerFactory =  vtkSmartPointer<vtkImageReader2Factory>::New();
                vtkImageReader2 * imageReader = readerFactory->CreateImageReader2(path_to_image_file.c_str());
                imageReader->SetFileName(path_to_image_file.c_str());
                imageReader->Update();
                vtk_image = imageReader->GetOutput();
            }
        }

        int extent[6];
        vtk_image->GetExtent(extent);
        double spacing[3];
        vtk_image->GetSpacing(spacing);
        double origin[3];
        vtk_image->GetOrigin(origin);


        LOG(DEBUG) << "VTK-Image information: ";
        LOG(DEBUG) << "-- Origin:  " << origin[0] << ", " << origin[1] << ", " << origin[2];
        LOG(DEBUG) << "-- Spacing: " << spacing[0] << ", " << spacing[1] << ", " << spacing[2];
        LOG(DEBUG) << "-- Extent:  " << "x = [" << extent[0] << ", " << extent[1] << "], y = [ " << extent[2] << ", " << extent[3] << "], z = [ " << extent[4] << ", " << extent[5] << "]";

    }
    else
    {
        LOG(FATAL) << "File '" << path_to_image_file << "' does not exist.";
    }

    return vtk_image;
}


vtkSmartPointer<vtkImageData> loadVtkImageData(std::string path_to_file)
{
    vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
    // check file type
    std::string vti_suffix = ".vti";
    if (hf::hasSuffix(path_to_file, vti_suffix))
    {
        LOG(DEBUG)<< "Attempting to load image file '" << path_to_file <<"' as vtkImageData";
        vtkSmartPointer<vtkXMLImageDataReader> imageReader = vtkSmartPointer<vtkXMLImageDataReader>::New();
        imageReader->SetFileName(path_to_file.c_str());
        imageReader->Update();
        image = imageReader->GetOutput();
        LOG(DEBUG)<< "Loaded image file.";
    }
    else
    {
        LOG(FATAL) << "File '" << path_to_file <<"' does not have the expected file ending '.vti'.";
    }
    return image;
}

void writeVtkUnstructuredGrid(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, std::string path_to_output_file)
{
    LOG(DEBUG) << "Writing UnstructuredGrid to " << path_to_output_file;
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =   vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(path_to_output_file.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(unstructuredGrid);
#else
    writer->SetInputData(unstructuredGrid);
#endif
    writer->Write();
}


vtkSmartPointer<vtkUnstructuredGrid> loadVtkUnstructuredGrid(std::string path_to_mesh_file)
{
    vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    // check file type
    std::string vtu_suffix = ".vtu";
    if (hf::hasSuffix(path_to_mesh_file, vtu_suffix))
    {
        LOG(DEBUG)<< "Attempting to load mesh '" << path_to_mesh_file <<"' as unstructuredGrid with vtkXMLUnstructuredGridReader";
        vtkSmartPointer<vtkXMLUnstructuredGridReader> ugridReader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        ugridReader->SetFileName(path_to_mesh_file.c_str());
        ugridReader->Update();
        ugrid = ugridReader->GetOutput();
        LOG(DEBUG)<< "Loaded mesh";
    }
    else
    {
        LOG(FATAL) << "File '" << path_to_mesh_file <<"' does not have the expected file ending '.vtu'.";
    }
    return ugrid;
}




vtkSmartPointer<vtkPolyData> extractSurfaceFromVtkDataSet(vtkSmartPointer<vtkDataSet> vtkdataset)
{
    vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =  vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceFilter->SetInputData(vtkdataset);
    surfaceFilter->Update();
    vtkSmartPointer<vtkPolyData> vtkpolydata = surfaceFilter->GetOutput();
    LOG(DEBUG) << "Submesh Polydata has " << vtkpolydata->GetNumberOfPoints() << " points." ;
    LOG(DEBUG) << "Submesh Polydata has " << vtkpolydata->GetNumberOfVerts() << " vertices." ;
    LOG(DEBUG) << "Submesh Polydata has " << vtkpolydata->GetNumberOfPolys() << " polygons." ;
    LOG(DEBUG) << "Submesh Polydata has " << vtkpolydata->GetNumberOfCells() << " cells." ;

    return vtkpolydata;
}


std::list<int> getBlockIdList(vtkSmartPointer<vtkDataSet> dataset)
{
    //! Identify distinct values in "ElementBlockIds" Cell data array
    std::list<int> blockIdList = getUniqeIntList(dataset, "ElementBlockIds");

    return blockIdList;
}


std::list<int> getUniqeIntList(vtkSmartPointer<vtkDataSet> dataset, std::string array_name)
{
    //! Identify distinct values in "ElementBlockIds" Cell data array
    vtkSmartPointer<vtkCellData> cellData = dataset->GetCellData();
    vtkSmartPointer<vtkIntArray> intDataArray = vtkIntArray::SafeDownCast(cellData->GetArray(array_name.c_str()));
    std::list<int> int_list;
    for(vtkIdType i = 0; i < intDataArray->GetNumberOfTuples(); i++)
    {
        int_list.push_back(intDataArray->GetValue(i));
    }
    int_list.sort();
    int_list.unique();

    return int_list;
}

    std::list<vtkSmartPointer<vtkPolyData> > vtkPolyDataListFromVtkUnstructuredGrid(vtkSmartPointer<vtkUnstructuredGrid> vtk_unstructured_grid)
{
    std::list<int>  blockIdList = getBlockIdList(vtk_unstructured_grid);

    int index = 0;
    std::list<vtkSmartPointer<vtkPolyData> > vtkPolyData_list;
    for (std::list<int>::iterator it=blockIdList.begin(); it!=blockIdList.end(); ++it)
    {
        index ++;
        int blockId = *it;
        /// Split unstructured Grid into sub-grids
        LOG(DEBUG) << "splitting vtkunstructuredGrid into subgrids";
        LOG(DEBUG) << "Creating sub-mesh " << index << " for grey value " << blockId;
        vtkSmartPointer<vtkThreshold> selector = vtkSmartPointer<vtkThreshold>::New();
        selector->SetInputData(vtk_unstructured_grid);
        selector->SetInputArrayToProcess(0, 0, 0,
                                         vtkDataObject::FIELD_ASSOCIATION_CELLS,
                                         "ElementBlockIds");
        selector->ThresholdBetween(blockId,blockId);
        selector->Update();
        vtkSmartPointer<vtkUnstructuredGrid> subgrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        subgrid = selector->GetOutput();
        LOG(DEBUG) << "Number of cells in subgrid: " << subgrid->GetNumberOfCells() << std::endl;

        /// extract surfaces from subgrid
        LOG(DEBUG) << "extracting surface from subgrid" << index;
        vtkSmartPointer<vtkPolyData> vtkpolydata = extractSurfaceFromVtkDataSet(subgrid);
        LOG(DEBUG) << "extracting surface from subgrid" << index;

        /// Assign BlockId as CellData
        LOG(DEBUG) << "BlockId " << blockId << " as CellData";
        vtkSmartPointer<vtkCellData> cellData = vtkpolydata->GetCellData();
        cellData->SetActiveScalars("ElementBlockIds");
        double * range = cellData->GetScalars()->GetRange();
        int labelValue = (int)(*range);
        std::cout<< "Grey label value = " << labelValue << std::endl;
        vtkPolyData_list.push_back(vtkpolydata);
    }

    return vtkPolyData_list;
}


vtkSmartPointer<vtkImageData> createVtkImageFromVtkPolyDataList(std::list<vtkSmartPointer<vtkPolyData> > vtkPolyDataList, double bounds[6], double spacing[3])
{
    vtkSmartPointer<vtkImageWeightedSum> sumFilter = vtkSmartPointer<vtkImageWeightedSum>::New();
    sumFilter->NormalizeByWeightOff();

    LOG(DEBUG)<< "Create vtkImage for each vtkPolyData" ;
    int index = 0;
    for (std::list<vtkSmartPointer<vtkPolyData> >::iterator it=vtkPolyDataList.begin(); it!=vtkPolyDataList.end(); ++it)
    {
        vtkSmartPointer<vtkPolyData> polydata = *it;
        vtkSmartPointer<vtkImageData> vtkImage = createVtkImageFromVtkPolyData(polydata, bounds, spacing);

        LOG(DEBUG) << "Adding image " << index << " to ImageWeightedSum filter";
        sumFilter->SetWeight(index,1);
        sumFilter->AddInputData(vtkImage);

//         std::string path_to_image_file = generatePathToFile("", "segment", ".nii", index);
//         writeVtkImage(vtkImage, path_to_image_file);
//         std::string path_to_output_file = generatePathToFile("", "segment", ".vtp", index);
//         writeVtkPolyData(polydata, path_to_output_file);
        index ++;
    }

    sumFilter->Update();
    vtkSmartPointer<vtkImageCast> summedCastFilter =  vtkSmartPointer<vtkImageCast>::New();
    summedCastFilter->SetInputConnection(sumFilter->GetOutputPort());
    summedCastFilter->SetOutputScalarTypeToUnsignedChar();
    summedCastFilter->Update();

    vtkSmartPointer<vtkImageData> vtk_image = summedCastFilter->GetOutput();
    return vtk_image;
}



vtkSmartPointer<vtkImageData> createVtkImageFromVtkPolyData(vtkSmartPointer<vtkPolyData> vtkPolyData, double bounds[6], double spacing[3])
{
    LOG(DEBUG) << "Getting BlockId";
    vtkSmartPointer<vtkCellData> cellData = vtkPolyData->GetCellData();
    cellData->SetActiveScalars("ElementBlockIds");
    double * range = cellData->GetScalars()->GetRange();
    int labelValue = (int)(*range);
    LOG(DEBUG)<< "Grey label value = " << labelValue ;

    // create Image
    vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();
    whiteImage->SetSpacing(spacing);

    LOG(DEBUG) << "Start image creation";
    // compute dimensions
    int dim[3];
    for (int i = 0; i < 3; i++)
    {
        dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]));
    }
    whiteImage->SetDimensions(dim);
    //whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);
    whiteImage->SetExtent(0, dim[0], 0, dim[1], 0, dim[2]);
    double origin[3];
    origin[0] = bounds[0] + spacing[0] / 2;
    origin[1] = bounds[2] + spacing[1] / 2;
    origin[2] = bounds[4] + spacing[2] / 2;
    whiteImage->SetOrigin(origin);

#if VTK_MAJOR_VERSION <= 5
    whiteImage->SetScalarTypeToUnsignedChar();
whiteImage->AllocateScalars();
#else
    whiteImage->AllocateScalars(VTK_UNSIGNED_INT,1);
#endif
    // fill the image with foreground voxels for current blockID:
    int inval = labelValue;
    LOG(DEBUG) << " Assigned pixel foreground value :" << inval << std::endl;
    vtkIdType count = whiteImage->GetNumberOfPoints();
    LOG(DEBUG) << "number of points in image " << count;
    for (vtkIdType i = 0; i < count; ++i)
    {
        whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
    }

    // polygonal data --> image stencil:
    vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc =   vtkSmartPointer<vtkPolyDataToImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
    pol2stenc->SetInput(polydata);
#else
    pol2stenc->SetInputData(vtkPolyData);
#endif
    pol2stenc->SetOutputOrigin(origin);
    pol2stenc->SetOutputSpacing(spacing);
    pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
    pol2stenc->Update();

    // cut the corresponding white image and set the background:
    unsigned char outval = 0;
    vtkSmartPointer<vtkImageStencil> imgstenc =  vtkSmartPointer<vtkImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
    imgstenc->SetInput(whiteImage);
#else
    imgstenc->SetInputData(whiteImage);
    imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
#endif
    imgstenc->ReverseStencilOff();
    imgstenc->SetBackgroundValue(outval);
    imgstenc->Update();

    vtkSmartPointer<vtkImageData> vtkImage = imgstenc->GetOutput();

    return vtkImage;
}


vtkSmartPointer<vtkImageData> createVtkImageFromVtkUnstructuredGrid(vtkSmartPointer<vtkUnstructuredGrid> vtk_ugrid, double spacing[3])
{
    double input_bounds[6];
    vtk_ugrid->GetBounds(input_bounds);
    LOG(DEBUG) << "Bounds of input mesh: " << input_bounds[0] << ", "
               << input_bounds[1] << ", "
               << input_bounds[2] << ", "
               << input_bounds[3] << ", "
               << input_bounds[4] << ", "
               << input_bounds[5] << ".";

    std::list<vtkSmartPointer<vtkPolyData> > vtkPolyData_list = vtkPolyDataListFromVtkUnstructuredGrid(vtk_ugrid);
    vtkSmartPointer<vtkImageData> vtk_image = createVtkImageFromVtkPolyDataList(vtkPolyData_list, input_bounds, spacing);
    return vtk_image;
}



void writeVtkImageData(vtkSmartPointer<vtkImageData> image_data, std::string path_to_output_file)
{
    LOG(DEBUG) << "Writing imageData to " << path_to_output_file;
    vtkSmartPointer<vtkXMLImageDataWriter> writer =   vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(path_to_output_file.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(image_data);
#else
    writer->SetInputData(image_data);
#endif
    writer->Write();
}


void writeVtkData(vtkSmartPointer<vtkDataSet> dataset, std::string path_to_output_file)
{   std::string vtp_suffix = ".vtp";
    std::string vti_suffix = ".vti";
    std::string vtu_suffix = ".vtu";

    if (hf::hasSuffix(path_to_output_file, vtp_suffix))
    {
        LOG(DEBUG) << "Writing vtkPolyData to '" << path_to_output_file << "'";
        vtkSmartPointer<vtkPolyData> data = vtkPolyData::SafeDownCast(dataset);
        vtkSmartPointer<vtkXMLPolyDataWriter> writer =   vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName(path_to_output_file.c_str());
        writer->SetInputData(data);
        writer->Write();
    }
    else if (hf::hasSuffix(path_to_output_file, vti_suffix))
    {
        LOG(DEBUG) << "Writing vtkImageData to '" << path_to_output_file << "'";
        vtkSmartPointer<vtkImageData> data = vtkImageData::SafeDownCast(dataset);
        vtkSmartPointer<vtkXMLImageDataWriter> writer =   vtkSmartPointer<vtkXMLImageDataWriter>::New();
        writer->SetFileName(path_to_output_file.c_str());
        writer->SetInputData(data);
        writer->Write();
    }
    else if (hf::hasSuffix(path_to_output_file, vtu_suffix))
    {
        LOG(DEBUG) << "Writing vtkUnstructuredGrid to '" << path_to_output_file << "'";
        vtkSmartPointer<vtkUnstructuredGrid> data = vtkUnstructuredGrid::SafeDownCast(dataset);
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =   vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(path_to_output_file.c_str());
        writer->SetInputData(data);
        writer->Write();
    }
    else
    {
        LOG(DEBUG) << "Cannot write '" << path_to_output_file << "' -- Unknown file ending.";
    }

}

}  // end helperftcs
