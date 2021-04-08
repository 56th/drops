/// \author Alexander Zhiliakov alex@math.uh.edu

#ifndef VTK_WRITER_HPP
#define VTK_WRITER_HPP

#ifdef _VTK
    #include "vtkSmartPointer.h"
    #include "vtkPoints.h"
    #include "vtkLagrangeTetra.h"
    #include "vtkCellArray.h"
    #include "vtkUnstructuredGrid.h"
    #include "vtkXMLUnstructuredGridWriter.h"
    #include "vtkDoubleArray.h"
    #include "vtkPointData.h"
#endif
#include "geom/multigrid.h"
#include "SingletonLogger.hpp"

namespace DROPS {
    class VTKWriter {
    public:
        struct VTKVar {
            std::string name;
            VecDescCL& vec;
        };
        using GridBuilder = void(VTKWriter::*)();
        struct VTKParams {
            std::string path;
            MultiGridCL const * mg = nullptr;
            bool binary = true;
            bool staticGrid = true;
            GridBuilder builder = &VTKWriter::buildFullGrid;
        };
    private:
        VTKParams const & params;
        std::vector<VTKVar> vars;
        std::unordered_map<VertexCL const *, size_t> vertexIndex;
        std::unordered_map<EdgeCL const *, size_t> edgeIndex;
        size_t frame = 0;
        #ifdef _VTK
            vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
        #endif
    public:
        VTKWriter(VTKParams const & params) : params(params) {}
        void buildFullGrid() { // update mesh, cf. https://lorensen.github.io/VTKExamples/site/Cxx/IO/WriteVTU/
            #ifdef _VTK
                vertexIndex.clear();
                edgeIndex.clear();
                unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
                size_t n = 0;
                auto points = vtkSmartPointer<vtkPoints>::New();
                for (auto it = params.mg->GetTriangVertexBegin(); it != params.mg->GetTriangVertexEnd(); ++it) {
                    points->InsertNextPoint(it->GetCoord()[0], it->GetCoord()[1], it->GetCoord()[2]);
                    vertexIndex[&*it] = n++;
                }
                for (auto it = params.mg->GetTriangEdgeBegin(); it != params.mg->GetTriangEdgeEnd(); ++it) {
                    auto p = GetBaryCenter(*it);
                    points->InsertNextPoint(p[0], p[1], p[2]);
                    edgeIndex[&*it] = n++;
                }
                auto cellArray = vtkSmartPointer<vtkCellArray>::New();
                for (auto it = params.mg->GetTriangTetraBegin(); it != params.mg->GetTriangTetraEnd(); ++it) {
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
                unstructuredGrid->SetPoints(points);
                unstructuredGrid->SetCells(VTK_LAGRANGE_TETRAHEDRON, cellArray);
            #endif
        }
        void buildInterfaceGrid() {
            #ifdef _VTK
                vertexIndex.clear();
                edgeIndex.clear();
                unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
                auto points = vtkSmartPointer<vtkPoints>::New();
                auto cellArray = vtkSmartPointer<vtkCellArray>::New();
                size_t n = 0;
                auto isActive = [&](TetraCL const & tet) {
                    for (auto const & var : vars)
                        if (var.vec.RowIdx->IsOnInterface() && tet.Unknowns.Exist(var.vec.RowIdx->GetIdx()))
                            return true;
                    return false;
                };
                for (auto it = params.mg->GetTriangTetraBegin(); it != params.mg->GetTriangTetraEnd(); ++it) {
                    if (!isActive(*it)) continue;
                    for (auto vit = it->GetVertBegin(); vit != it->GetVertEnd(); ++vit)
                        if (!vertexIndex.count(*vit)) {
                            points->InsertNextPoint((*vit)->GetCoord()[0], (*vit)->GetCoord()[1], (*vit)->GetCoord()[2]);
                            vertexIndex[*vit] = n++;
                        }
                    for (auto eit = it->GetEdgesBegin(); eit != it->GetEdgesEnd(); ++eit)
                        if (!edgeIndex.count(*eit)) {
                            auto p = GetBaryCenter(**eit);
                            points->InsertNextPoint(p[0], p[1], p[2]);
                            edgeIndex[*eit] = n++;
                        }
                    auto tetra = vtkSmartPointer<vtkLagrangeTetra>::New();
                    tetra->GetPointIds()->SetNumberOfIds(10);
                    tetra->GetPoints()->SetNumberOfPoints(10);
                    tetra->Initialize();
                    for (size_t i : {0, 1, 2, 3})
                        tetra->GetPointIds()->SetId(i, vertexIndex[it->GetVertex(i)]);
                    tetra->GetPointIds()->SetId(4, edgeIndex[it->GetEdge(EdgeByVert(0, 1))]);
                    tetra->GetPointIds()->SetId(5, edgeIndex[it->GetEdge(EdgeByVert(1, 2))]);
                    tetra->GetPointIds()->SetId(6, edgeIndex[it->GetEdge(EdgeByVert(0, 2))]);
                    tetra->GetPointIds()->SetId(7, edgeIndex[it->GetEdge(EdgeByVert(0, 3))]);
                    tetra->GetPointIds()->SetId(8, edgeIndex[it->GetEdge(EdgeByVert(1, 3))]);
                    tetra->GetPointIds()->SetId(9, edgeIndex[it->GetEdge(EdgeByVert(2, 3))]);
                    cellArray->InsertNextCell(tetra);
                }
                unstructuredGrid->SetPoints(points);
                unstructuredGrid->SetCells(VTK_LAGRANGE_TETRAHEDRON, cellArray);
            #endif
        }
        size_t numVars() const { return vars.size(); }
        VTKWriter& add(VTKVar const & var) {
            vars.push_back(var);
            return *this;
        }
        VTKWriter& write(double time) {
            std::string funcName = __func__;
            if (!params.mg) throw std::invalid_argument(funcName + ": provide mg");
            #ifdef _VTK
                auto& logger = SingletonLogger::instance();
                if (!unstructuredGrid || !params.staticGrid) {
                    logger.beg("build grid");
                        (this->*params.builder)();
                    logger.end();
                }
                logger.beg("prepare vars");
                    auto n = unstructuredGrid->GetNumberOfPoints();
                    for (auto& var : vars) {
                        logger.beg("var: " + var.name);
                            size_t dim = var.vec.RowIdx->IsScalar() ? 1 : 3;
                            if (var.vec.RowIdx->NumUnknowns() % dim) throw std::logic_error(funcName + ": invalid numb of d.o.f.");
                            auto size = dim * n;
                            auto blockSize = var.vec.RowIdx->NumUnknowns() / dim;
                            auto idx = var.vec.RowIdx->GetIdx();
                            auto array = vtkSmartPointer<vtkDoubleArray>::New();
                            logger.buf
                                << "FE space idx: " << idx << '\n'
                                << "dim: " << dim << '\n'
                                << "size: " << size;
                            logger.log();
                            auto* value = new double[size];
                            for (size_t i = 0; i < size; ++i) value[i] = 0.;
                            for (auto it = params.mg->GetTriangVertexBegin(); it != params.mg->GetTriangVertexEnd(); ++it) {
                                auto i = &*it;
                                if (vertexIndex.count(i) && it->Unknowns.Exist(idx))
                                    for (size_t d = 0; d < dim; ++d)
                                        value[dim * vertexIndex[i] + d] = var.vec.Data[it->Unknowns(idx) + d * blockSize];
                            }
                            for (auto it = params.mg->GetTriangEdgeBegin(); it != params.mg->GetTriangEdgeEnd(); ++it) {
                                auto i = &*it;
                                if (edgeIndex.count(i)) {
                                    if (it->Unknowns.Exist(idx))
                                        for (size_t d = 0; d < dim; ++d)
                                            value[dim * edgeIndex[i] + d] = var.vec.Data[it->Unknowns(idx) + d * blockSize];
                                    else { // P1
                                        auto i0 = dim * vertexIndex[it->GetVertex(0)];
                                        auto i1 = dim * vertexIndex[it->GetVertex(1)];
                                        for (size_t d = 0; d < dim; ++d)
                                            value[dim * edgeIndex[i] + d] = .5 * (value[i0 + d] + value[i1 + d]);
                                    }
                                }
                            }
                            array->SetNumberOfComponents(dim); // cf. https://vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataCellNormals
                            array->SetArray(value, size, 0 /* this will free the memory for ptr, cf. https://vtk.org/doc/release/5.6/html/a00505.html#d35ae5bee4aa873d543f1ab3eaf94454 */);
                            array->SetName(var.name.c_str());
                            unstructuredGrid->GetPointData()->AddArray(array);
                        logger.end();
                    }
                logger.end();
                logger.beg("write .vtu");
                    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
                    auto name = params.path + '_' + std::to_string(frame) + ".vtu";
                    writer->SetFileName(name.c_str());
                    writer->SetInputData(unstructuredGrid);
                    if (!params.binary) writer->SetDataModeToAscii();
                    writer->Write();
                logger.end();
                logger.beg("write .pvd");
                    name = name.substr(name.find_last_of("/\\") + 1);
                    if  (frame == 0) {
                        std::ofstream pvd(params.path + ".pvd");
                        pvd <<
                            "<?xml version=\"1.0\"?>\n"
                            "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
                            "<Collection>\n"
                                "\t<DataSet timestep=\""<< time <<"\" group=\"\" part=\"0\" file=\"" << name <<"\"/>\n"
                            "</Collection>\n"
                            "</VTKFile>";
                    } else {
                        std::ofstream pvd;
                        pvd.open(params.path + ".pvd", std::ios_base::in);
                        if(!pvd.is_open()) std::logic_error(funcName + ": cannot reopen .pvd file");
                        pvd.seekp(-24, std::ios_base::end);
                        pvd <<
                                "\t<DataSet timestep=\""<< time <<"\" group=\"\" part=\"0\" file=\"" << name <<"\"/>\n"
                            "</Collection>\n"
                            "</VTKFile>";
                        pvd.close();
                    }
                logger.end();
                ++frame;
                return *this;
            #else
                throw std::logic_error(funcName + ": DROPS compiled w/o VTK lib");
            #endif
        }
    };
}

#endif // VTK_WRITER_HPP