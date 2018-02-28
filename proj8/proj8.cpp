#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkDataSetWriter.h>

int main(){
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj8.vtk");

    //Contour the data
    vtkCountourFilter *cf = vtkCountourFilter::New();
    cf->SetNumberOfContours(1);
    cf->SetValue(0, 3.0);
    cf->SetInputConnection(rdr->GetOutputPort());

    vtkDataSetWriter *wrtr = vtkDataSetWriter::New();
    wrtr->SetFileName("contour.vtk");
    wrtr->SetInputConnection(cf->GetOutputPort());
    wrtr->Write();

    rdr->Delete();
    cf->Delete();
    wrter->Delete();
    
}