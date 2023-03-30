/*

Copyright (c) 2005-2018, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "DomMeshBasedCellPopulationWithGhostNodes.hpp"
#include "CellLocationIndexWriter.hpp"
#include "VtkMeshWriter.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeVelocityWriter.hpp"
#include "CellPopulationElementWriter.hpp"
#include "CellId.hpp"
#include "Debug.hpp"
#include "MutableMesh.hpp"

#include "WildTypeCellMutationState.hpp"
#include "StromalCellMutationState.hpp"


template<unsigned DIM>
DomMeshBasedCellPopulationWithGhostNodes<DIM>::DomMeshBasedCellPopulationWithGhostNodes(
     MutableMesh<DIM, DIM>& rMesh,
     std::vector<CellPtr>& rCells,
     const std::vector<unsigned> locationIndices,
     bool deleteMesh,
     double ghostSpringStiffness,
     double ghostRestSeperation)
             : MeshBasedCellPopulation<DIM,DIM>(rMesh, rCells, locationIndices, deleteMesh, false), // do not call the base class Validate()
            //    mpVoronoiTessellation(nullptr),
               mGhostSpringStiffness(ghostSpringStiffness),
               mGhostRestSeperation(ghostRestSeperation),
               mWriteVtkAsPointsDom(false),
               mOutputMeshInVtkDom(false)
{

    if (!locationIndices.empty())
    {
        // Create a set of node indices corresponding to ghost nodes
        std::set<unsigned> node_indices;
        std::set<unsigned> location_indices;
        std::set<unsigned> ghost_node_indices;

        for (unsigned i=0; i<this->GetNumNodes(); i++)
        {
            node_indices.insert(this->GetNode(i)->GetIndex());
        }
        for (unsigned i=0; i<locationIndices.size(); i++)
        {
            location_indices.insert(locationIndices[i]);
        }

        std::set_difference(node_indices.begin(), node_indices.end(),
                            location_indices.begin(), location_indices.end(),
                            std::inserter(ghost_node_indices, ghost_node_indices.begin()));

        // This method finishes and then calls Validate()
        SetGhostNodes(ghost_node_indices);
    }
    else
    {
        this->mIsGhostNode = std::vector<bool>(this->GetNumNodes(), false);
        Validate();
    }

}

template<unsigned DIM>
DomMeshBasedCellPopulationWithGhostNodes<DIM>::DomMeshBasedCellPopulationWithGhostNodes(MutableMesh<DIM, DIM>& rMesh,
                                                                                  double ghostSpringStiffness, double ghostRestSeperation)
    : MeshBasedCellPopulation<DIM,DIM>(rMesh),
      mGhostSpringStiffness(ghostSpringStiffness),
      mGhostRestSeperation(ghostRestSeperation)
{
}

template<unsigned DIM>
DomMeshBasedCellPopulationWithGhostNodes<DIM>::~DomMeshBasedCellPopulationWithGhostNodes()
{
}

template<unsigned DIM>
TetrahedralMesh<DIM, DIM>* DomMeshBasedCellPopulationWithGhostNodes<DIM>::GetTetrahedralMeshForPdeModifier()
{
    EXCEPTION("Currently can't solve PDEs on meshes with ghost nodes");
    return static_cast<TetrahedralMesh<DIM, DIM>*>(&(this->mrMesh));
}

template<unsigned DIM>
std::vector<bool>& DomMeshBasedCellPopulationWithGhostNodes<DIM>::rGetGhostNodes()
{
    return this->mIsGhostNode;
}

template<unsigned DIM>
bool DomMeshBasedCellPopulationWithGhostNodes<DIM>::IsGhostNode(unsigned index)
{
    return this->mIsGhostNode[index];
}

template<unsigned DIM>
std::set<unsigned> DomMeshBasedCellPopulationWithGhostNodes<DIM>::GetGhostNodeIndices()
{
    std::set<unsigned> ghost_node_indices;
    for (unsigned i=0; i<this->mIsGhostNode.size(); i++)
    {
        if (this->mIsGhostNode[i])
        {
            ghost_node_indices.insert(i);
        }
    }
    return ghost_node_indices;
}

template<unsigned DIM>
void DomMeshBasedCellPopulationWithGhostNodes<DIM>::SetGhostNodes(const std::set<unsigned>& rGhostNodeIndices)
{
    // Reinitialise all entries of mIsGhostNode to false
    this->mIsGhostNode = std::vector<bool>(this->mrMesh.GetNumNodes(), false);

    // Update mIsGhostNode
    for (std::set<unsigned>::iterator iter=rGhostNodeIndices.begin(); iter!=rGhostNodeIndices.end(); ++iter)
    {
        this->mIsGhostNode[*iter] = true;
    }

    Validate();
}

template<unsigned DIM>
c_vector<double, DIM> DomMeshBasedCellPopulationWithGhostNodes<DIM>::CalculateForceBetweenGhostNodes(const unsigned& rNodeAGlobalIndex, const unsigned& rNodeBGlobalIndex)
{
    assert(rNodeAGlobalIndex != rNodeBGlobalIndex);
    c_vector<double, DIM> unit_difference;
    c_vector<double, DIM> force_between_nodes;
    const c_vector<double, DIM>& r_node_a_location = this->GetNode(rNodeAGlobalIndex)->rGetLocation();
    const c_vector<double, DIM>& r_node_b_location = this->GetNode(rNodeBGlobalIndex)->rGetLocation();

    // There is reason not to subtract one position from the other (cylindrical meshes)
    unit_difference = this->mrMesh.GetVectorFromAtoB(r_node_a_location, r_node_b_location);

    double distance_between_nodes = norm_2(unit_difference);
    unit_difference /= distance_between_nodes;

    double rest_length;
    if(this->IsGhostNode(rNodeAGlobalIndex) && this->IsGhostNode(rNodeBGlobalIndex))
    {
        // rest_length = mGhostRestSeperation;
        rest_length = 1.0;
    }
    else
    {
        // rest_length = 1.0;
        rest_length = mGhostRestSeperation;
    }
    
    double cut_off = 1.5;
    if (distance_between_nodes <= cut_off)
    {
        force_between_nodes = mGhostSpringStiffness * unit_difference * (distance_between_nodes - rest_length);
    }
    else
    {
        force_between_nodes[0] = 0.0;
        force_between_nodes[1] = 0.0;
        force_between_nodes[2] = 0.0;
    }
    return force_between_nodes;
}

template<unsigned DIM>
CellPtr DomMeshBasedCellPopulationWithGhostNodes<DIM>::AddCell(CellPtr pNewCell, CellPtr pParentCell)
{
    // Add new cell to population
    CellPtr p_created_cell = MeshBasedCellPopulation<DIM,DIM>::AddCell(pNewCell, pParentCell);
    assert(p_created_cell == pNewCell);

    // Update size of mIsGhostNode if necessary
    unsigned new_node_index = this->GetLocationIndexUsingCell(p_created_cell);

    if (this->GetNumNodes() > this->mIsGhostNode.size())
    {
        this->mIsGhostNode.resize(this->GetNumNodes());
        this->mIsGhostNode[new_node_index] = false;
    }

    // Return pointer to new cell
    return p_created_cell;
}

template<unsigned DIM>
void DomMeshBasedCellPopulationWithGhostNodes<DIM>::Validate()
{
    // Get a list of all the nodes that are ghosts
    std::vector<bool> validated_node = mIsGhostNode;
    assert(mIsGhostNode.size()==this->GetNumNodes());

    // Look through all of the cells and record what node they are associated with.
    for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = this->GetLocationIndexUsingCell((*cell_iter));

        // If the node attached to this cell is labelled as a ghost node, then throw an error
        if (mIsGhostNode[node_index])
        {
            EXCEPTION("Node " << node_index << " is labelled as a ghost node and has a cell attached");
        }
        validated_node[node_index] = true;
    }

    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if (!validated_node[i])
        {
            EXCEPTION("Node " << i << " does not appear to be a ghost node or have a cell associated with it");
        }
    }
}

template<unsigned DIM>
void DomMeshBasedCellPopulationWithGhostNodes<DIM>::UpdateGhostNodesAfterReMesh(NodeMap& rMap)
{
    
    // Copy mIsGhostNode to a temporary vector
    std::vector<bool> ghost_nodes_before_remesh = mIsGhostNode;

    // Reinitialise mIsGhostNode
    mIsGhostNode.clear();
    mIsGhostNode.resize(this->GetNumNodes());

    // Update mIsGhostNode using the node map
    for (unsigned old_index=0; old_index<rMap.GetSize(); old_index++)
    {
        if (!rMap.IsDeleted(old_index))
        {
            unsigned new_index = rMap.GetNewIndex(old_index);
            mIsGhostNode[new_index] = ghost_nodes_before_remesh[old_index];
        }
    }
}

template<unsigned DIM>
std::set<unsigned> DomMeshBasedCellPopulationWithGhostNodes<DIM>::GetNeighbouringLocationIndices(CellPtr pCell)
{
    unsigned node_index = this->GetLocationIndexUsingCell(pCell);
    std::set<unsigned> neighbour_indices = this->GetNeighbouringNodeIndices(node_index);

    // Remove ghost nodes from the neighbour indices
    for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
         iter != neighbour_indices.end();)
    {
        if (this->IsGhostNode(*iter))
        {
            neighbour_indices.erase(iter++);
        }
        else
        {
            ++iter;
        }
    }

    return neighbour_indices;
}

template<unsigned DIM>
void DomMeshBasedCellPopulationWithGhostNodes<DIM>::AcceptCellWritersAcrossPopulation()
{
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = this->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        // If it isn't a ghost node then there might be cell writers attached
        if (! this->IsGhostNode(node_iter->GetIndex()))
        {
            for (typename std::vector<boost::shared_ptr<AbstractCellWriter<DIM, DIM> > >::iterator cell_writer_iter = this->mCellWriters.begin();
                 cell_writer_iter != this->mCellWriters.end();
                 ++cell_writer_iter)
            {
                CellPtr cell_from_node = this->GetCellUsingLocationIndex(node_iter->GetIndex());
                this->AcceptCellWriter(*cell_writer_iter, cell_from_node);
            }
        }
    }
}

template<unsigned DIM>
void DomMeshBasedCellPopulationWithGhostNodes<DIM>::ApplyGhostForces()
{
    // Initialise vector of forces on ghost nodes
    std::vector<c_vector<double, DIM> > drdt(this->GetNumNodes());
    for (unsigned i=0; i<drdt.size(); i++)
    {
        drdt[i] = zero_vector<double>(DIM);
    }

    // Calculate forces on ghost nodes
    for (typename MutableMesh<DIM, DIM>::EdgeIterator edge_iterator = static_cast<MutableMesh<DIM, DIM>&>((this->mrMesh)).EdgesBegin();
        edge_iterator != static_cast<MutableMesh<DIM, DIM>&>((this->mrMesh)).EdgesEnd();
        ++edge_iterator)
    {
        unsigned nodeA_global_index = edge_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = edge_iterator.GetNodeB()->GetIndex();

        c_vector<double, DIM> force = CalculateForceBetweenGhostNodes(nodeA_global_index, nodeB_global_index);

        if (!this->mIsGhostNode[nodeA_global_index])
        {
            if (this->mIsGhostNode[nodeB_global_index])
            {
                drdt[nodeB_global_index] -= force;
            }
        }
        else if (this->mIsGhostNode[nodeA_global_index])
        {
            drdt[nodeA_global_index] += force;

            if (this->mIsGhostNode[nodeB_global_index])
            {
                drdt[nodeB_global_index] -= force;
            }
        }
    }

    for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
         node_iter != this->mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        if (this->mIsGhostNode[node_index])
        {
            node_iter->ClearAppliedForce();
            node_iter->AddAppliedForceContribution(drdt[node_index]);
        }
    }
}

template<unsigned DIM>
void DomMeshBasedCellPopulationWithGhostNodes<DIM>::OpenWritersFiles(OutputFileHandler& rOutputFileHandler)
{
    if (this->mOutputResultsForChasteVisualizer)
    {
        if (!this-> template HasWriter<CellLocationIndexWriter>())
        {
            this-> template AddCellWriter<CellLocationIndexWriter>();
        }
    }

    MeshBasedCellPopulation<DIM, DIM>::OpenWritersFiles(rOutputFileHandler);
}

template<unsigned DIM>
void DomMeshBasedCellPopulationWithGhostNodes<DIM>::WriteVtkResultsToFile(const std::string& rDirectory)
{
#ifdef CHASTE_VTK
// PRINT_VARIABLE((this->mpVoronoiTessellation != nullptr));
    if (this->mpVoronoiTessellation != nullptr)
    {
        unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
        std::stringstream time;
        time << num_timesteps;

        // Create mesh writer for VTK output
        VertexMeshWriter<DIM, DIM> mesh_writer(rDirectory, "results", false);

        // Iterate over any cell writers that are present
        unsigned num_vtk_cells = this->mpVoronoiTessellation->GetNumElements();
        for (typename std::vector<boost::shared_ptr<AbstractCellWriter<DIM, DIM> > >::iterator cell_writer_iter = this->mCellWriters.begin();
             cell_writer_iter != this->mCellWriters.end();
             ++cell_writer_iter)
        {
            // Create vector to store VTK cell data
            std::vector<double> vtk_cell_data(num_vtk_cells);

            // Loop over elements of mpVoronoiTessellation
            for (typename VertexMesh<DIM, DIM>::VertexElementIterator elem_iter = this->mpVoronoiTessellation->GetElementIteratorBegin();
                 elem_iter != this->mpVoronoiTessellation->GetElementIteratorEnd();
                 ++elem_iter)
            {
                // Get the indices of this element and the corresponding node in mrMesh
                unsigned elem_index = elem_iter->GetIndex();
                unsigned node_index = this->mpVoronoiTessellation->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);

                // If this node corresponds to a ghost node, set any "cell" data to be -1.0
                if (this->IsGhostNode(node_index))
                {
                    // Populate the vector of VTK cell data
                    vtk_cell_data[elem_index] = -1.0;
                }
                else
                {
                    // Get the cell corresponding to this node
                    CellPtr p_cell = this->GetCellUsingLocationIndex(node_index);

                    // Populate the vector of VTK cell data
                    vtk_cell_data[elem_index] = (*cell_writer_iter)->GetCellDataForVtkOutput(p_cell, this);
                }
            }

            mesh_writer.AddCellData((*cell_writer_iter)->GetVtkCellDataName(), vtk_cell_data);
        }

        // Next, record which nodes are ghost nodes
        // Note that the cell writer hierarchy can not be used to do this as ghost nodes don't have corresponding cells.
        std::vector<double> ghosts(num_vtk_cells);
        for (typename VertexMesh<DIM, DIM>::VertexElementIterator elem_iter = this->mpVoronoiTessellation->GetElementIteratorBegin();
             elem_iter != this->mpVoronoiTessellation->GetElementIteratorEnd();
             ++elem_iter)
        {
            unsigned elem_index = elem_iter->GetIndex();
            unsigned node_index = this->mpVoronoiTessellation->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);
            ghosts[elem_index]  = (double) (this->IsGhostNode(node_index));
        }
        mesh_writer.AddCellData("Non-ghosts", ghosts);

        ///\todo #1975 - deal with possibility of information stored in CellData

        mesh_writer.WriteVtkUsingMesh(*(this->mpVoronoiTessellation), time.str());
        *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
        *(this->mpVtkMetaFile) << num_timesteps;
        *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
        *(this->mpVtkMetaFile) << num_timesteps;
        *(this->mpVtkMetaFile) << ".vtu\"/>\n";
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Store the present time as a string
    unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    std::stringstream time;
    time << num_timesteps;

    // When outputting any CellData, we assume that the first cell is representative of all cells
    unsigned num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
    std::vector<std::string> cell_data_names = this->Begin()->GetCellData()->GetKeys();

    //PRINT_VECTOR(cell_data_names);
    // std::vector<std::vector<double> > cell_data;

    // if (mOutputMeshInVtkDom)
    // {
    //     // Create mesh writer for VTK output
    //     // VtkMeshWriter<DIM, DIM> mesh_writer(rDirectory, "mesh_"+time.str(), false);
    //     // mesh_writer.WriteFilesUsingMesh(rGetMesh());
    //     // PRINT_VARIABLE(rGetMesh().GetNumNodes());
    //     // PRINT_VARIABLE(rGetMesh().GetNumElements());

    //     // Creat a new mesh
    //     MutableMesh<2,3>* pMutableMesh = new MutableMesh<2,3>();

    //     // Store the nodes from the old mesh
    //     std::vector<Node<3>*> nodes;

    //     double dist_threshold = 2.0;

    //     for (unsigned node_index=0; node_index<rGetMesh().GetNumNodes(); node_index++)
    //     {
    //         const c_vector<double, 3>& r_location = rGetMesh().GetNode(node_index)->rGetLocation();
    //         nodes.push_back(new Node<3>(node_index, r_location));
    //         pMutableMesh->AddNode(new Node<3>(node_index, r_location));
    //     }


    //     int element_iter = 0;
    //     for (unsigned elem_index=0; elem_index<rGetMesh().GetNumElements(); elem_index++)
    //     {    

    //         Element<DIM,DIM>* pElement=rGetMesh().GetElement(elem_index);

    //         // PRINT_VARIABLE(pElement);

    //         for( int j=0; j<1; j++) //for( int j=0; j<4; j++)
    //         {
    //             unsigned Node0Index = pElement->GetNodeGlobalIndex((0+j)%4);
    //             unsigned Node1Index = pElement->GetNodeGlobalIndex((1+j)%4);
    //             unsigned Node2Index = pElement->GetNodeGlobalIndex((2+j)%4);
    //             unsigned Node3Index = pElement->GetNodeGlobalIndex((3+j)%4);

    //             const c_vector<double, 3>& r_location0 = rGetMesh().GetNode(Node0Index)->rGetLocation();
    //             const c_vector<double, 3>& r_location1 = rGetMesh().GetNode(Node1Index)->rGetLocation();
    //             const c_vector<double, 3>& r_location2 = rGetMesh().GetNode(Node2Index)->rGetLocation();
    //             const c_vector<double, 3>& r_location3 = rGetMesh().GetNode(Node3Index)->rGetLocation();


    //             double distance_between_nodes_01 = norm_2(r_location0 - r_location1);
    //             double distance_between_nodes_02 = norm_2(r_location0 - r_location2);
    //             double distance_between_nodes_03 = norm_2(r_location0 - r_location3);
    //             double distance_between_nodes_12 = norm_2(r_location1 - r_location2);
    //             double distance_between_nodes_13 = norm_2(r_location1 - r_location3);
    //             double distance_between_nodes_23 = norm_2(r_location2 - r_location3);

    //             if( distance_between_nodes_01 <= dist_threshold && distance_between_nodes_03 <= dist_threshold && distance_between_nodes_13 <= dist_threshold )
    //             {
    //                 std::vector<Node<3>*> ElementNodes0;
    //                 ElementNodes0.push_back(nodes[Node0Index]);
    //                 ElementNodes0.push_back(nodes[Node1Index]);
    //                 ElementNodes0.push_back(nodes[Node3Index]);
    //                 pMutableMesh->AddElement(new Element<2,3>(element_iter,ElementNodes0));
    //                 element_iter++;
    //             }

    //             if( distance_between_nodes_01 <= dist_threshold && distance_between_nodes_02 <= dist_threshold && distance_between_nodes_12 <= dist_threshold )
    //             {
    //                 std::vector<Node<3>*> ElementNodes1;
    //                 ElementNodes1.push_back(nodes[Node0Index]);
    //                 ElementNodes1.push_back(nodes[Node1Index]);
    //                 ElementNodes1.push_back(nodes[Node2Index]);
    //                 pMutableMesh->AddElement(new Element<2,3>(element_iter,ElementNodes1));
    //                 element_iter++;
    //             }

    //             if( distance_between_nodes_02 <= dist_threshold && distance_between_nodes_03 <= dist_threshold && distance_between_nodes_23 <= dist_threshold )
    //             {
    //                 std::vector<Node<3>*> ElementNodes2;
    //                 ElementNodes2.push_back(nodes[Node0Index]);
    //                 ElementNodes2.push_back(nodes[Node2Index]);
    //                 ElementNodes2.push_back(nodes[Node3Index]);
    //                 pMutableMesh->AddElement(new Element<2,3>(element_iter,ElementNodes2));
    //                 element_iter++;
    //             }

    //             if( distance_between_nodes_12 <= dist_threshold && distance_between_nodes_13 <= dist_threshold && distance_between_nodes_23 <= dist_threshold )
    //             {
    //                 std::vector<Node<3>*> ElementNodes3;
    //                 ElementNodes3.push_back(nodes[Node1Index]);
    //                 ElementNodes3.push_back(nodes[Node2Index]);
    //                 ElementNodes3.push_back(nodes[Node3Index]);
    //                 pMutableMesh->AddElement(new Element<2,3>(element_iter,ElementNodes3));
    //                 element_iter++;
    //             }

    //         }
    //     }
        
    //     std::ostringstream time_string;
    //     time_string << SimulationTime::Instance()->GetTimeStepsElapsed();
    //     std::string results_file = "mesh_results_" + time_string.str();
    //     VtkMeshWriter<2,3>* p_vtk_mesh_writer = new VtkMeshWriter<2,3>(rDirectory, results_file, false);
    //     p_vtk_mesh_writer->WriteFilesUsingMesh(*pMutableMesh);

    //     delete pMutableMesh;
    // }

    if (mWriteVtkAsPointsDom)
    {
        //std::cout<< "\n";
        //PRINT_VARIABLE("Start ");
        
        // Create mesh writer for VTK output
        VtkMeshWriter<DIM, DIM> cells_writer(rDirectory, "results_"+time.str(), false);

        
        // Iterate over any cell writers that are present
        unsigned num_cells = this->GetNumAllCells();
        unsigned total_num_nodes = this->mrMesh.GetNumNodes();
        //PRINT_2_VARIABLES(total_num_nodes, num_cell_data_items);
        
        unsigned cell_data[num_cell_data_items][total_num_nodes];

        std::vector<double> ghosts(total_num_nodes);

        for (typename std::vector<boost::shared_ptr<AbstractCellWriter<DIM, DIM> > >::iterator cell_writer_iter = this->mCellWriters.begin();
             cell_writer_iter != this->mCellWriters.end();
             ++cell_writer_iter)
        {
            // Create vector to store VTK cell data
            std::vector<double> vtk_cell_data(total_num_nodes);

            // Loop over cells
            for (unsigned node_index=0; node_index<total_num_nodes; node_index++)
            {
                // If this node corresponds to a ghost node, set any "cell" data to be -1.0
                if (this->IsGhostNode(node_index))
                {
                    // Populate the vector of VTK cell data
                    vtk_cell_data[node_index] = -1.0;
                }
                else
                {
                    // Get the cell corresponding to this node
                    CellPtr p_cell = this->GetCellUsingLocationIndex(node_index);

                    // Populate the vector of VTK cell data
                    vtk_cell_data[node_index] = (*cell_writer_iter)->GetCellDataForVtkOutput(p_cell, this);
                    //PRINT_VARIABLE((*cell_writer_iter)->GetCellDataForVtkOutput(p_cell, this));
                }
                
            }

            cells_writer.AddPointData((*cell_writer_iter)->GetVtkCellDataName(), vtk_cell_data);
        }
        // Loop over cells
        for (unsigned node_index=0; node_index<total_num_nodes; node_index++)
        {

            // If this node corresponds to a ghost node, set any "cell" data to be -1.0
            if (this->IsGhostNode(node_index))
            {
                for (unsigned var=0; var<num_cell_data_items; var++)
                {
                    cell_data[var][node_index] = -1;
                    
                }
            }
            else
            {
                // Get the cell corresponding to this node
                CellPtr p_cell = this->GetCellUsingLocationIndex(node_index);

                for (unsigned var=0; var<num_cell_data_items; var++)
                {
                    // TRACE("adding data");
                    cell_data[var][node_index] = p_cell->GetCellData()->GetItem(cell_data_names[var]);
                    // TRACE("got data");
                }
            }
            
        }
        
        // for (unsigned var=0; var<num_cell_data_items; var++)
        // {
        //     cells_writer.AddPointData(cell_data_names[var], cell_data[var]);
        // }
       
        // Make a copy of the nodes in a disposable mesh for writing
        {
            std::vector<Node<DIM>* > nodes;
            for (unsigned index=0; index<this->mrMesh.GetNumNodes(); index++)
            {
            //
            //PRINT_VARIABLE((this->IsGhostNode(index)));
            //
            Node<DIM>* p_node = this->mrMesh.GetNode(index);
            nodes.push_back(p_node);
            
            /*    if(!(this->IsGhostNode(index)))
                {
                    //TRACE(this->mrMesh.GetNode(index));
                    Node<DIM>* p_node = this->mrMesh.GetNode(index);
                    nodes.push_back(p_node);
                }
                else if (this->IsGhostNode(index))
                {
                    //Do something...
                }
                //PRINT_2_VARIABLES(index,(this->IsGhostNode(index)));
            */
            }
            NodesOnlyMesh<DIM> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5); // Arbitrary cut off as connectivity not used.
            cells_writer.WriteFilesUsingMesh(mesh);
        }
        //PRINT_2_VARIABLES("Middle",num_timesteps);
        
        
        *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
        *(this->mpVtkMetaFile) << num_timesteps;
        *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
        *(this->mpVtkMetaFile) << num_timesteps;
        *(this->mpVtkMetaFile) << ".vtu\"/>\n";
        
        //PRINT_VARIABLE("Finish");
    }

    // if (false)//Change this to true to display nodes without ghosts, but it needs fixing....
    // {
    //     //std::cout<< "\n";
    //     //PRINT_VARIABLE("Start ");
        
    //     // Create mesh writer for VTK output
    //     VtkMeshWriter<DIM, DIM> cells_writer(rDirectory, "results_no_ghost"+time.str(), false);

        
    //     // Iterate over any cell writers that are present
    //     unsigned num_cells = this->GetNumAllCells();
    //     unsigned total_num_nodes = this->mrMesh.GetNumNodes();
    //     //PRINT_2_VARIABLES(total_num_nodes, num_cell_data_items);
        
    //     std::vector<double> ghosts(total_num_nodes);

    //     for (typename std::vector<boost::shared_ptr<AbstractCellWriter<DIM, DIM> > >::iterator cell_writer_iter = this->mCellWriters.begin();
    //          cell_writer_iter != this->mCellWriters.end();
    //          ++cell_writer_iter)
    //     {
    //         // Create vector to store VTK cell data
    //         std::vector<double> vtk_cell_data(total_num_nodes);

    //         // Loop over cells
    //         for (unsigned node_index=0; node_index<total_num_nodes; node_index++)
    //         {
    //             // If this node corresponds to a ghost node, set any "cell" data to be -1.0
    //             if (this->IsGhostNode(node_index))
    //             {
    //                 // Populate the vector of VTK cell data
    //                 //vtk_cell_data[node_index] = -1.0;
    //             }
    //             else
    //             {
    //                 // Get the cell corresponding to this node
    //                 CellPtr p_cell = this->GetCellUsingLocationIndex(node_index);

    //                 // Populate the vector of VTK cell data
    //                 vtk_cell_data[node_index] = (*cell_writer_iter)->GetCellDataForVtkOutput(p_cell, this);
    //                 //PRINT_VARIABLE((*cell_writer_iter)->GetCellDataForVtkOutput(p_cell, this));

    //                 cells_writer.AddPointData((*cell_writer_iter)->GetVtkCellDataName(), vtk_cell_data);
    //             }
                
    //         }

            
    //     }
        
    //     // Loop over cells
    //     for (unsigned node_index=0; node_index<total_num_nodes; node_index++)
    //     {

    //         // If this node corresponds to a ghost node, set any "cell" data to be -1.0
    //         if (this->IsGhostNode(node_index))
    //         {
    //             for (unsigned var=0; var<num_cell_data_items; var++)
    //             {
    //                 //cell_data[var][node_index] = -1;
                    
    //             }
    //         }
    //         else
    //         {
    //             // Get the cell corresponding to this node
    //             CellPtr p_cell = this->GetCellUsingLocationIndex(node_index);

    //             for (unsigned var=0; var<num_cell_data_items; var++)
    //             {
    //                 cell_data[var][node_index] = p_cell->GetCellData()->GetItem(cell_data_names[var]);

    //                 cells_writer.AddPointData(cell_data_names[var], cell_data[var]);
                    
    //             }
    //         }
            
    //     }
        
    //     // for (unsigned var=0; var<num_cell_data_items; var++)
    //     // {
    //     //     cells_writer.AddPointData(cell_data_names[var], cell_data[var]);
    //     // }
       
    //     // Make a copy of the nodes in a disposable mesh for writing
    //     {
    //         std::vector<Node<DIM>* > nodes;
    //         for (unsigned index=0; index<this->mrMesh.GetNumNodes(); index++)
    //         {
            
    //         // PRINT_VARIABLE((this->IsGhostNode(index)));
            
    //         if ( (this->IsGhostNode(index)) )
    //         {
    //             Node<DIM>* p_node = this->mrMesh.GetNode(index);

    //             nodes.push_back(p_node);
    //         }
    //         else
    //         {
    //             Node<DIM>* p_node = this->mrMesh.GetNode(index);
    //             nodes.push_back(p_node);
    //         }
    //         /*    if(!(this->IsGhostNode(index)))
    //             {
    //                 //TRACE(this->mrMesh.GetNode(index));
    //                 Node<DIM>* p_node = this->mrMesh.GetNode(index);
    //                 nodes.push_back(p_node);
    //             }
    //             else if (this->IsGhostNode(index))
    //             {
    //                 //Do something...
    //             }
    //             //PRINT_2_VARIABLES(index,(this->IsGhostNode(index)));
    //         */
    //         }
    //         NodesOnlyMesh<DIM> mesh;
    //         mesh.ConstructNodesWithoutMesh(nodes, 1.5); // Arbitrary cut off as connectivity not used.
    //         cells_writer.WriteFilesUsingMesh(mesh);
    //     }
    //     //PRINT_2_VARIABLES("Middle",num_timesteps);
        
        
    //     *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    //     *(this->mpVtkMetaFile) << num_timesteps;
    //     *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    //     *(this->mpVtkMetaFile) << num_timesteps;
    //     *(this->mpVtkMetaFile) << ".vtu\"/>\n";
        
    //     //PRINT_VARIABLE("Finish");
    // }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif //CHASTE_VTK
}

////
template<unsigned DIM>
void DomMeshBasedCellPopulationWithGhostNodes<DIM>::SetWriteVtkAsPointsDom(bool writeVtkAsPointsDom)
{
    mWriteVtkAsPointsDom = writeVtkAsPointsDom;
}

template<unsigned DIM>
bool DomMeshBasedCellPopulationWithGhostNodes<DIM>::GetWriteVtkAsPointsDom()
{
    return mWriteVtkAsPointsDom;
}

template<unsigned DIM>
void DomMeshBasedCellPopulationWithGhostNodes<DIM>::SetOutputMeshInVtkDom(bool outputMeshInVtkDom)
{
    mOutputMeshInVtkDom = outputMeshInVtkDom;
}

template<unsigned DIM>
bool DomMeshBasedCellPopulationWithGhostNodes<DIM>::GetOutputMeshInVtkDom()
{
    return mOutputMeshInVtkDom;
}

/**
 * The cylindrical mesh is only defined in 2D, hence there is
 * a separate definition for this method in 3D, which doesn't have the capability
 * of dealing with periodic boundaries in 3D. This is \todo #1374.
 */
// template<>
// void DomMeshBasedCellPopulationWithGhostNodes<3>::CreateVoronoiTessellation()
// {
//     delete mpVoronoiTessellation;
//     mpVoronoiTessellation = new VertexMesh<3, 3>(static_cast<MutableMesh<3, 3> &>((this->mrMesh)));
// }

// template<unsigned DIM>
// VertexMesh<DIM,DIM>* DomMeshBasedCellPopulationWithGhostNodes<DIM>::GetVoronoiTessellation()
// {
//     assert(mpVoronoiTessellation!=nullptr);
//     return mpVoronoiTessellation;
// }

////
template<unsigned DIM>
void DomMeshBasedCellPopulationWithGhostNodes<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<GhostSpringStiffness>" << mGhostSpringStiffness << "</GhostSpringStiffness>\n";
    *rParamsFile << "\t\t<GhostRestSeperation>" << mGhostRestSeperation << "</GhostRestSeperation>\n";
    *rParamsFile << "\t\t<WriteVtkAsPoints>" << mWriteVtkAsPointsDom << "</WriteVtkAsPoints>\n";
    *rParamsFile << "\t\t<OutputMeshInVtk>" << mOutputMeshInVtkDom << "</OutputMeshInVtk>\n";

    // Call method on direct parent class
    MeshBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

// Explicit instantiation
template class DomMeshBasedCellPopulationWithGhostNodes<1>;
template class DomMeshBasedCellPopulationWithGhostNodes<2>;
template class DomMeshBasedCellPopulationWithGhostNodes<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DomMeshBasedCellPopulationWithGhostNodes)
