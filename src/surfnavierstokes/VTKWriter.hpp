/// \author Alexander Zhiliakov alex@math.uh.edu

#ifndef VTK_WRITER_HPP
#define VTK_WRITER_HPP

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkQuadraticTetra.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
//#include <vtkVertexGlyphFilter.h>

#include "geom/multigrid.h"

namespace DROPS {
    class VTKWriter {
    private:
        struct VTKVar {
            std::string name;
            VectorCL* value;
        };
        std::vector<VTKVar> vars;
        MultiGridCL const * mg;
        std::string path;
        std::unordered_map<VertexCL const *, size_t> vertexIndex;
        std::unordered_map<EdgeCL const *, size_t> edgeIndex;
    public:
        VTKWriter(std::string const & path, MultiGridCL const & mg) : path(path), mg(&mg) {}
        VTKWriter& add(std::string const & name, VectorCL& value) {
            vars.push_back({ name, &value });
            return *this;
        }
        VTKWriter& write(double t) {
            // (1) update mesh, cf. https://lorensen.github.io/VTKExamples/site/Cxx/IO/WriteVTU/
            vertexIndex.clear();
            edgeIndex.clear();
            size_t n = 0;
            auto points = vtkSmartPointer<vtkPoints>::New();
            for (auto it = mg->GetTriangVertexBegin(); it != mg->GetTriangVertexEnd(); ++it) {
                points->InsertNextPoint(it->GetCoord()[0], it->GetCoord()[1], it->GetCoord()[2]);
                vertexIndex[&*it] = n++;
            }
            std::cout << "numb of vertices = " << n << '\n';
            for (auto it = mg->GetTriangEdgeBegin(); it != mg->GetTriangEdgeEnd(); ++it) {
                auto p = GetBaryCenter(*it);
                points->InsertNextPoint(p[0], p[1], p[2]);
                edgeIndex[&*it] = n++;
            }
            std::cout << "numb of vertices + edges = " << n << '\n';
            auto cellArray = vtkSmartPointer<vtkCellArray>::New();
            for (auto it = mg->GetTetrasBegin(); it != mg->GetTetrasEnd(); ++it) {
                auto tetra = vtkSmartPointer<vtkQuadraticTetra>::New();
                for (size_t i : {0, 1, 2, 3})
                    tetra->GetPointIds()->SetId(i, vertexIndex[it->GetVertex(i)]);
                // cf. https://vtk.org/doc/nightly/html/classvtkQuadraticTetra.html#details
                tetra->GetPointIds()->SetId(4, edgeIndex[it->GetEdge(EdgeByVert(0, 1))]);
                tetra->GetPointIds()->SetId(5, edgeIndex[it->GetEdge(EdgeByVert(1, 2))]);
                tetra->GetPointIds()->SetId(6, edgeIndex[it->GetEdge(EdgeByVert(0, 2))]);
                tetra->GetPointIds()->SetId(7, edgeIndex[it->GetEdge(EdgeByVert(0, 3))]);
                tetra->GetPointIds()->SetId(8, edgeIndex[it->GetEdge(EdgeByVert(1, 3))]);
                tetra->GetPointIds()->SetId(9, edgeIndex[it->GetEdge(EdgeByVert(2, 3))]);
                cellArray->InsertNextCell(tetra);
            }
            auto unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
            unstructuredGrid->SetPoints(points);
            unstructuredGrid->SetCells(VTK_QUADRATIC_TETRA, cellArray);
            // (2) update vars
            for (auto const & var : vars) {
                auto array = vtkSmartPointer<vtkDoubleArray>::New();
                array->SetArray(&var.value->operator[](0), var.value->size(), 1);
                std::cout << "var size = " << var.value->size() << '\n';
                array->SetName(var.name.c_str());
                unstructuredGrid->GetPointData()->AddArray(array);
            }
            // (3) write
            auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
            writer->SetFileName((path + ".vtu").c_str());
            writer->SetInputData(unstructuredGrid);
            writer->Write();
            return *this;
        }
    };
}

#endif // VTK_WRITER_HPP