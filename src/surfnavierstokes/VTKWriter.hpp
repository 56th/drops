/// \author Alexander Zhiliakov alex@math.uh.edu

#ifndef VTK_WRITER_HPP
#define VTK_WRITER_HPP

#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkLagrangeTetra.h"
#include "vtkCellArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"

#include "geom/multigrid.h"

namespace DROPS {
    class VTKWriter {
    public:
        struct VTKVar {
            enum class Type { P2, vecP2, P1, vecP1 };
            std::string name;
            VectorCL* value;
            Type type;
        };
    private:
        bool binary;
        std::vector<VTKVar> vars;
        MultiGridCL const * mg;
        std::string path;
        std::unordered_map<VertexCL const *, size_t> vertexIndex;
        std::unordered_map<EdgeCL const *, size_t> edgeIndex;
        size_t frame = 0;
    public:
        VTKWriter(std::string const & path, MultiGridCL const & mg, bool binary = true) : path(path), mg(&mg), binary(binary) {}
        VTKWriter& add(VTKVar const & var) {
            vars.push_back(var);
            return *this;
        }
        VTKWriter& write(double t) {
            std::string funcName = __func__;
            // (1) update mesh, cf. https://lorensen.github.io/VTKExamples/site/Cxx/IO/WriteVTU/
            vertexIndex.clear();
            edgeIndex.clear();
            size_t n = 0;
            auto points = vtkSmartPointer<vtkPoints>::New();
            for (auto it = mg->GetTriangVertexBegin(); it != mg->GetTriangVertexEnd(); ++it) {
                points->InsertNextPoint(it->GetCoord()[0], it->GetCoord()[1], it->GetCoord()[2]);
                vertexIndex[&*it] = n++;
            }
            for (auto it = mg->GetTriangEdgeBegin(); it != mg->GetTriangEdgeEnd(); ++it) {
                auto p = GetBaryCenter(*it);
                points->InsertNextPoint(p[0], p[1], p[2]);
                edgeIndex[&*it] = n++;
            }
            auto cellArray = vtkSmartPointer<vtkCellArray>::New();
            for (auto it = mg->GetTetrasBegin(); it != mg->GetTetrasEnd(); ++it) {
                auto tetra = vtkSmartPointer<vtkLagrangeTetra>::New();
                tetra->GetPointIds()->SetNumberOfIds(10);
                tetra->GetPoints()->SetNumberOfPoints(10);
                tetra->Initialize();
                for (size_t i : {0, 1, 2, 3})
                    tetra->GetPointIds()->SetId(i, vertexIndex[it->GetVertex(i)]);
                // for nodes ordering, cf. https://vtk.org/doc/nightly/html/classvtkQuadraticTetra.html#details and https://blog.kitware.com/modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/
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
            unstructuredGrid->SetCells(VTK_LAGRANGE_TETRAHEDRON, cellArray);
            // (2) update vars
            for (auto const & var : vars) {
                if(var.type == VTKVar::Type::P2) {
                    if (var.value->size() != n)
                        throw std::invalid_argument(funcName + ": inconsistent number of d.o.f. for P2 interpolant");
                    auto array = vtkSmartPointer<vtkDoubleArray>::New();
                    array->SetArray(&var.value->operator[](0), n, 1);
                    array->SetName(var.name.c_str());
                    unstructuredGrid->GetPointData()->AddArray(array);
                }
                if (var.type == VTKVar::Type::P1) {
                    auto  array = vtkSmartPointer<vtkDoubleArray>::New();
                    std::vector<double> value;
                    value.assign(&var.value->operator[](0), &var.value->operator[](0) + var.value->size());
                    for (auto it = mg->GetTriangEdgeBegin(); it != mg->GetTriangEdgeEnd(); ++it) {
                        auto i1 = vertexIndex[it->GetVertex(0)];
                        auto i2 = vertexIndex[it->GetVertex(1)];
                        value.push_back(.5 * (value[i1] + value[i2]));
                    }
                    if (value.size() != n) throw std::invalid_argument(funcName + ": inconsistent number of d.o.f. for P1 interpolant");
                    auto* ptr = new double[n];
                    std::copy(value.begin(), value.end(), ptr);
                    array->SetArray(ptr, n, 0 /* this will free the memory for ptr, cf. https://vtk.org/doc/release/5.6/html/a00505.html#d35ae5bee4aa873d543f1ab3eaf94454 */);
                    array->SetName(var.name.c_str());
                    unstructuredGrid->GetPointData()->AddArray(array);
                }
                if (var.type == VTKVar::Type::vecP2) {
                    if (var.value->size() != 3 * n)
                        throw std::invalid_argument(funcName + ": inconsistent number of d.o.f. for vecP2 interpolant");
                    auto array = vtkSmartPointer<vtkDoubleArray>::New();
                    array->SetNumberOfComponents(3); // cf. https://vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataCellNormals
                    array->SetArray(&var.value->operator[](0), 3 * n, 1);
                    array->SetName(var.name.c_str());
                    unstructuredGrid->GetPointData()->AddArray(array);
                }
                if (var.type == VTKVar::Type::vecP1)
                    throw std::logic_error(funcName + ": vecP1 export not yet implemented");
            }
            // (3) write
            auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
            writer->SetFileName((path + '_' + std::to_string(frame++) + ".vtu").c_str());
            writer->SetInputData(unstructuredGrid);
            if (!binary) writer->SetDataModeToAscii();
            writer->Write();
            return *this;
        }
    };
}

#endif // VTK_WRITER_HPP