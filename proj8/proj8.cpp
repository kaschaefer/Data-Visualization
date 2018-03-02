#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkDataSetWriter.h>
#include <vtkPolyDataMapper.h>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkDataSetMapper.h>
#include <vtkCamera.h>
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkProperty.h>
#include <vtkHedgeHog.h>
#include <vtkPointData.h>
#include <vtkDataSetAttributes.h>
#include <vtkThreshold.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkDataSet.h>
#include <vtkExtractRectilinearGrid.h>
#include <vtkStreamTracer.h>
#include <vtkPlaneSource.h>
#include <vtkLineSource.h>


#include <iostream>

int main(int argc, char* argv[]){
    //Get DataSet Reader
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj8.vtk");
    rdr->Update();

    //Get HardyGlobal Data
    vtkDataSet *hardyGlobal = rdr->GetOutput();
    hardyGlobal->GetPointData()->SetActiveScalars("hardyglobal");

    //***Data Flow for Renderer 1
    //Configure Contour Filter
    vtkContourFilter *cf = vtkContourFilter::New();
    cf->SetNumberOfContours(2);
    cf->SetValue(0, 2.5);
    cf->SetValue(1, 5.0);
    cf->SetInputData(hardyGlobal);
    cf->Update();

    //Map Data to Filter
    vtkSmartPointer<vtkDataSetMapper> contourMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    contourMapper->SetInputConnection(cf->GetOutputPort());
    contourMapper->SetScalarRange(rdr->GetOutput()->GetPointData()->GetScalars()->GetRange());

    //***Data Flow for Renderer 2
    //Configure Complete 3D Data Set
    vtkSmartPointer<vtkDataSetMapper> rectMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    rectMapper->SetInputData(hardyGlobal);
    rectMapper->Update();
    
    //Create Planes
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    vtkSmartPointer<vtkPlane> plane2 = vtkSmartPointer<vtkPlane>::New();
    vtkSmartPointer<vtkPlane> plane3 = vtkSmartPointer<vtkPlane>::New();
    
    plane->SetOrigin(0,0,0);
    plane->SetNormal(0,0,1);
    plane2->SetOrigin(0,0,0);
    plane2->SetNormal(0,1,0);
    plane3->SetOrigin(0,0,0);
    plane3->SetNormal(1,0,0);

    //Configure Cut Filters to Cut Planes
    vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
    vtkSmartPointer<vtkCutter> cutter2 = vtkSmartPointer<vtkCutter>::New();
    vtkSmartPointer<vtkCutter> cutter3 = vtkSmartPointer<vtkCutter>::New();

    cutter->SetCutFunction(plane);
    cutter->SetInputData(hardyGlobal);
    cutter2->SetCutFunction(plane2);
    cutter2->SetInputData(hardyGlobal);
    cutter3->SetCutFunction(plane3);
    cutter3->SetInputData(hardyGlobal);

    cutter->Update();
    cutter2->Update();
    cutter3->Update();

    //Map Data to Filter
    vtkSmartPointer<vtkDataSetMapper> cutterMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    vtkSmartPointer<vtkDataSetMapper> cutterMapper2 = vtkSmartPointer<vtkDataSetMapper>::New();
    vtkSmartPointer<vtkDataSetMapper> cutterMapper3 = vtkSmartPointer<vtkDataSetMapper>::New();    
    cutterMapper->SetInputConnection(cutter->GetOutputPort());
    cutterMapper2->SetInputConnection(cutter2->GetOutputPort());
    cutterMapper3->SetInputConnection(cutter3->GetOutputPort());
    cutterMapper->SetScalarRange(rdr->GetOutput()->GetPointData()->GetScalars()->GetRange());
    cutterMapper2->SetScalarRange(rdr->GetOutput()->GetPointData()->GetScalars()->GetRange());
    cutterMapper3->SetScalarRange(rdr->GetOutput()->GetPointData()->GetScalars()->GetRange());


    //grad is a field F
    //Hedgehog filter for Renderer 3
    //hedgehog filter -- needs to be pointed at the data
    //Set active attribute
    vtkDataSet *grad = rdr->GetOutput();
    grad->GetPointData()->SetActiveVectors("grad");
    vtkExtractRectilinearGrid *extrtRectGrid = vtkExtractRectilinearGrid::New();
    extrtRectGrid->SetInputData(grad);
    extrtRectGrid->SetSampleRate(5, 5, 5);
    extrtRectGrid->Update();
    //Create Filters
    vtkHedgeHog *hedgehog = vtkHedgeHog::New();
    vtkArrowSource *arrowSource = vtkArrowSource::New();
    vtkGlyph3D *glyph = vtkGlyph3D::New();
    //Configure HH Filter
    hedgehog->SetInputData(extrtRectGrid->GetOutput());
    hedgehog->SetScaleFactor(5);
    hedgehog->Update();
    //Configure ArrowSource
    arrowSource->Update();
    //Configure Glyphs
    glyph->SetInputData(hedgehog->GetOutput());
    glyph->SetSourceConnection(arrowSource->GetOutputPort());
    glyph->OrientOn();
    glyph->Update();

    vtkSmartPointer<vtkDataSetMapper> hedgeHogMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    hedgeHogMapper->SetInputConnection(glyph->GetOutputPort());
    hedgeHogMapper->SetScalarRange(grad->GetScalarRange());

    vtkSmartPointer<vtkActor> hedgehogActor = vtkSmartPointer<vtkActor>::New();
    hedgehogActor->SetMapper(hedgeHogMapper);

    //StreamTracer filter for Renderer 4
    vtkStreamTracer *tracer = vtkStreamTracer::New();
    tracer->SetIntegratorType(vtkStreamTracer::RUNGE_KUTTA4);
    rdr->GetOutput()->GetPointData()->SetActiveVectors("grad");
    tracer->SetInputData(rdr->GetOutput());
    vtkLineSource *lineSrc = vtkLineSource::New();
    lineSrc->SetPoint1(-9,0,0);
    lineSrc->SetPoint2(0,0,-9);
    
    vtkPoints *points = vtkPoints::New();
    points->InsertNextPoint(-9,0,0);
    points->InsertNextPoint(-8,0,0);
    points->InsertNextPoint(-7,0,0);
    points->InsertNextPoint(-6,0,0);
    points->InsertNextPoint(-5,0,0);
    points->InsertNextPoint(-4,0,0);
    points->InsertNextPoint(-3,0,0);
    points->InsertNextPoint(-2,0,0);
    points->InsertNextPoint(-1,0,0);
    points->InsertNextPoint(0,0,0);
    points->InsertNextPoint(1,0,0);
    points->InsertNextPoint(2,0,0);
    points->InsertNextPoint(3,0,0);
    points->InsertNextPoint(4,0,0);
    points->InsertNextPoint(5,0,0);
    points->InsertNextPoint(7,0,0);
    points->InsertNextPoint(8,0,0);
    points->InsertNextPoint(9,0,0);

    lineSrc->SetPoints(points);
    tracer->SetSourceConnection(lineSrc->GetOutputPort());
    tracer->SetMaximumPropagation(100);
    tracer->SetInitialIntegrationStep(0.1);
    tracer->Update();

    vtkSmartPointer<vtkDataSetMapper> tracerMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    tracerMapper->SetInputConnection(tracer->GetOutputPort());
    tracerMapper->Update();
    tracerMapper->SetScalarRange(grad->GetScalarRange());
    

    //Actors
        //Renderer 1
    vtkSmartPointer<vtkActor> ren1Actor = vtkSmartPointer<vtkActor>::New();
    ren1Actor->SetMapper(contourMapper);

        //Renderer 2
    vtkSmartPointer<vtkActor> planeActor = vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkActor> plane2Actor = vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkActor> plane3Actor = vtkSmartPointer<vtkActor>::New();

    planeActor->GetProperty()->SetColor(255,255,0);
    planeActor->GetProperty()->SetLineWidth(2);
    planeActor->SetMapper(cutterMapper);
    
    plane2Actor->GetProperty()->SetColor(255,255,0);
    plane2Actor->GetProperty()->SetLineWidth(2);
    plane2Actor->SetMapper(cutterMapper2);
    
    plane3Actor->GetProperty()->SetColor(255,255,0);
    plane3Actor->GetProperty()->SetLineWidth(2);
    plane3Actor->SetMapper(cutterMapper3);

    vtkSmartPointer<vtkActor> rectActor = vtkSmartPointer<vtkActor>::New();
    rectActor->GetProperty()->SetOpacity(0.5);
    rectActor->GetProperty()->SetColor(255,255,255);
    rectActor->SetMapper(rectMapper);

        //Renderer 4
    vtkSmartPointer<vtkActor> tracerActor = vtkSmartPointer<vtkActor>::New();
    tracerActor->SetMapper(tracerMapper);
    tracerActor->VisibilityOn();
        
    //Set Up Render Window and Interactor
    vtkSmartPointer<vtkRenderWindow> rendWindow = vtkSmartPointer<vtkRenderWindow>::New();
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();

    //Set Up Renderers
    vtkSmartPointer<vtkRenderer> ren1 = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderer> ren2 = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderer> ren3 = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderer> ren4 = vtkSmartPointer<vtkRenderer>::New();

    //Set Viewports
    ren1->SetViewport(0, 0, .5, .5);
    ren2->SetViewport(0, 0.5, 0.5, 1.0);
    ren3->SetViewport(0.5, 0, 1.0, 0.5);
    ren4->SetViewport(0.5, 0.5, 1.0, 1.0);

    //Add Renderers to Render Window
    rendWindow->AddRenderer(ren1);
    rendWindow->AddRenderer(ren2);
    rendWindow->AddRenderer(ren3);
    rendWindow->AddRenderer(ren4);

    //Add Interactor
    iren->SetRenderWindow(rendWindow);

    //Add Actors
    ren1->AddActor(ren1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);

    ren2->AddActor(planeActor);
    ren2->AddActor(plane2Actor);
    ren2->AddActor(plane3Actor);
    ren2->AddActor(rectActor);
    ren2->SetBackground(0,0,0);

    ren3->AddActor(hedgehogActor);
    ren3->SetBackground(0.0,0.0,0.0);

    ren4->AddActor(tracerActor);
    ren4->SetBackground(0.0,0.0,0.0);

    rendWindow->SetSize(600,600);

    //Configure Active Cameras
    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,70);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(70);

    ren2->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren2->GetActiveCamera()->SetPosition(0,0,70);
    ren2->GetActiveCamera()->SetViewUp(0,1,0);
    ren2->GetActiveCamera()->SetClippingRange(20, 120);
    ren2->GetActiveCamera()->SetDistance(70);

    ren3->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren3->GetActiveCamera()->SetPosition(0,0,70);
    ren3->GetActiveCamera()->SetViewUp(0,1,0);
    ren3->GetActiveCamera()->SetClippingRange(20, 120);
    ren3->GetActiveCamera()->SetDistance(70);

    //Invokes Initial Render
    iren->Initialize();
    iren->Start();

    return EXIT_SUCCESS;

}